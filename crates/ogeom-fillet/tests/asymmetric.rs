//! Asymmetric rounds: a different distance onto each face of a straight
//! edge between two planes, the section an ellipse or a conic.
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
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
    // The elliptical call reads the radius from the faces' legs, which may
    // stand an ulp off the closed form: the same solid to the bit at one of
    // the closed form's nearest floats.
    let want = (
        volume(&m, &asked.shape).to_bits(),
        vertex_bits(&m, &asked.shape),
    );
    let same = [radius.next_down(), radius, radius.next_up()]
        .into_iter()
        .any(|r| {
            let round = ogeom_fillet::fillet_edge(&mut m, &s, &edge, r, T).unwrap();
            (
                volume(&m, &round.shape).to_bits(),
                vertex_bits(&m, &round.shape),
            ) == want
        });
    assert!(same, "the circular round within an ulp of its radius");
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

    // Along a half drum's ruling, a leg longer than the flat is wide.
    let (mut m, d, side, flat, _) = half_drum();
    refused(
        e(
            &mut m,
            &d,
            &side,
            &flat,
            10.5,
            0.5,
            RoundProfile::Elliptical,
            T,
        ),
        "past that face's side",
    );
}

#[test]
fn blend_edges_takes_an_asymmetric_round_beside_others_and_refuses_a_mixed_pair_beside_a_concave_sharp_edge()
 {
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

    // A convex and a concave band at a vertex whose third edge, left
    // sharp, is concave: a block standing on the L's floor against its
    // wall, the round of its top's front edge and the fill of its top
    // against the wall, its upright against the wall left sharp. Refused
    // by name, asymmetric on either band.
    let mut m = Model::new();
    let (s, ..) = l_bracket(&mut m);
    let step = ogeom_algo::make_box(
        &mut m,
        Frame::new(Point::new(8.0, 0.0, 8.0), Direction::Z, Direction::X, T).unwrap(),
        (6.0, 6.0, 4.0),
        T,
    )
    .unwrap()
    .shape;
    let s = ogeom_bool::fuse(&mut m, &s, &step, T).unwrap().shape;
    let front = edge_where(&m, &s, |p| {
        (p.y - 6.0).abs() < 1e-9 && (p.z - 12.0).abs() < 1e-9
    });
    let back = edge_where(&m, &s, |p| {
        (p.x - 8.0).abs() < 1e-9 && (p.z - 12.0).abs() < 1e-9
    });
    let step_top = face_at(&m, &s, Point::new(11.0, 3.0, 12.0));
    for request in [
        [
            (front.clone(), asked(&step_top)),
            (back.clone(), Blend::Round(1.0)),
        ],
        [
            (front.clone(), Blend::Round(1.0)),
            (back.clone(), asked(&step_top)),
        ],
    ] {
        let message = ogeom_fillet::blend_edges(&mut m, &s, &request, ChamferCorner::Point, T)
            .expect_err("refused")
            .to_string();
        assert!(message.contains("both cut or both fill"), "{message}");
    }
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
fn level_laws_are_the_constant_asymmetric_round() {
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
}

/// A 10 high block on the quadrilateral (0, 0), (20, 0), (20 + lean, 20),
/// (0, 20): its end through (20, 0) leans back from the front edge by
/// `lean` over 20, or overhangs it where `lean` is negative.
fn leaning_block(m: &mut Model, lean: f64) -> Shape {
    let points = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(20.0, 0.0, 0.0),
        Point::new(20.0 + lean, 20.0, 0.0),
        Point::new(0.0, 20.0, 0.0),
    ];
    let wire = ogeom_algo::make_polygon(m, &points, true, T).unwrap().shape;
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

/// The top front edge of a block whose end at x = 20 leans, asked with
/// distances rising linearly along it, 1 to 2 on the top and 0.5 to 1.5 on
/// the front: the band's section at each station is the square conic
/// corner of that station's distances, so over the edge it takes the
/// integral of their product, `20 (1 - share) / 2 * int (1 + s)(0.5 + s)`.
/// Past the end the laws hold level. Where the end leans back the band
/// runs on until the end face cuts it, each point of its end section a
/// quarter of its depth into the top further: the end corner's area times
/// a quarter of its centroid's depth. Where the end overhangs it, the
/// band stops on it: at each station only the part of the section within
/// four times the station's distance from the end is under material.
#[test]
fn distances_varying_past_an_end_not_square_run_on_level_to_it() {
    use ogeom_fillet::RadiusLaw;
    let rho = 0.4;
    let w = conic(rho);
    let on_top = RadiusLaw::Linear {
        start: 1.0,
        end: 2.0,
    };
    let on_front = RadiusLaw::Linear {
        start: 0.5,
        end: 1.5,
    };
    // The share of its legs' triangle the square conic corner keeps, and
    // 20 int_0^t (1 + s)(0.5 + s) ds.
    let kept = corner_area(1.0, 1.0, FRAC_PI_2, w) / 0.5;
    let product = |t: f64| 20.0 * t * (0.5 + t * (0.75 + t / 3.0));
    let along = 0.5 * kept * product(1.0);
    for lean in [5.0, -5.0] {
        let mut m = Model::new();
        let s = leaning_block(&mut m, lean);
        let block = 10.0 * 20.0 * (40.0 + lean) * 0.5;
        let edge = edge_where(&m, &s, |p| p.y.abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9);
        let top = face_at(&m, &s, Point::new(10.0, 10.0, 10.0));
        let front = face_at(&m, &s, Point::new(10.0, 0.0, 5.0));
        let taken = if lean > 0.0 {
            let (area, centroid) = square_corner(2.0, 1.5, w);
            along + area * centroid[0] * lean / 20.0
        } else {
            // Up the top a distance y from the front, the end stands
            // -lean / 20 y short of x = 20: whole sections up to where the
            // end reaches the top's contact, then only what is under it,
            // integrated back from the end.
            let k = -lean / 20.0;
            let whole = (20.0 - k) / (1.0 + k / 20.0);
            0.5 * kept * product(whole / 20.0)
                + integrate(20.0 - whole, |back| {
                    let s = 1.0 - back / 20.0;
                    let (d_top, d_front) = (1.0 + s, 0.5 + s);
                    integrate(back / k, |y| conic_depth(d_top, d_front, w, y))
                })
        };
        let exact = block - taken;
        let mut volumes = Vec::new();
        for (face, a, b) in [(&top, &on_top, &on_front), (&front, &on_front, &on_top)] {
            let built = ogeom_fillet::fillet_edge_distances_variable(
                &mut m,
                &s,
                &edge,
                face,
                a,
                b,
                RoundProfile::Conic(rho),
                T,
            )
            .unwrap();
            let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{lean}: {:?}", diagnosis.problems);
            let v = volume(&m, &built.shape);
            assert!(
                (v - exact).abs() <= 1e-6,
                "{lean}: volume {v} against {exact}"
            );
            volumes.push(v);
        }
        assert_eq!(
            volumes[0].to_bits(),
            volumes[1].to_bits(),
            "naming the other face first"
        );
    }
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
        message.contains("left sharp") || message.contains("other than as one tangent chain"),
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

/// The box's corner at the origin's top: the top front edge along x, the
/// top left edge along y and the upright edge along z, with the faces the
/// requests name.
fn box_corner(m: &mut Model) -> (Shape, [Shape; 3], [Shape; 3]) {
    let (s, _, top, front) = block(m);
    let x = edge_where(m, &s, |p| p.y.abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9);
    let y = edge_where(m, &s, |p| p.x.abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9);
    let z = edge_where(m, &s, |p| p.x.abs() < 1e-9 && p.y.abs() < 1e-9);
    let left = face_at(m, &s, Point::new(0.0, 10.0, 5.0));
    (s, [x, y, z], [top, front, left])
}

/// The worst angle between a patch and the band faces it meets.
fn patch_tangency(m: &Model, s: &Shape) -> f64 {
    let mut worst = 0.0_f64;
    for f in explore_unique(m, s, ShapeType::Face).unwrap() {
        let fd = m.node(&f).unwrap().data().as_face().unwrap();
        if !matches!(
            m.geometry().surface(fd.surface),
            Some(ogeom_geom::SurfaceGeometry::BSpline(_))
        ) {
            continue;
        }
        for c in ogeom_fillet::analyse_blend(m, s, &f, 17, T).unwrap() {
            let nd = m.node(&c.neighbour).unwrap().data().as_face().unwrap();
            if matches!(
                m.geometry().surface(nd.surface),
                Some(
                    ogeom_geom::SurfaceGeometry::Extrusion(_)
                        | ogeom_geom::SurfaceGeometry::Cylinder(_)
                )
            ) {
                worst = worst.max(c.tangency_error);
            }
        }
    }
    worst
}

/// A box's convex corner where asymmetric rounds meet each other, a round
/// and a bevel: each band stopped at its setback and the opening closed by
/// one patch held inside the corner's planes, tangent to each round's
/// band. The same solid whichever face each asymmetric round names and in
/// whichever order the request lists its edges; and the corner takes no
/// more than its bands would along their whole edges, and more than the
/// widest band alone over its edge short of the corner.
#[test]
fn asymmetric_rounds_close_a_box_corner_with_rounds_and_bevels() {
    let asym = |face: &Shape, a: f64, b: f64| Blend::RoundDistances {
        face: face.clone(),
        on_face: a,
        on_other: b,
        profile: RoundProfile::Conic(0.45),
    };
    let w = conic(0.45);
    // Each request: the blends of x, y and z, given twice, naming the
    // other face the second time; and each band's section area.
    type Ask = fn(&[Shape; 3], bool, &dyn Fn(&Shape, f64, f64) -> Blend) -> [Blend; 3];
    let cases: [(&str, Ask, [f64; 3]); 3] = [
        (
            "three asymmetric rounds",
            |f, swap, asym| {
                if swap {
                    [
                        asym(&f[1], 1.0, 2.0),
                        asym(&f[2], 1.0, 1.5),
                        asym(&f[2], 0.8, 1.2),
                    ]
                } else {
                    [
                        asym(&f[0], 2.0, 1.0),
                        asym(&f[0], 1.5, 1.0),
                        asym(&f[1], 1.2, 0.8),
                    ]
                }
            },
            [
                corner_area(2.0, 1.0, FRAC_PI_2, conic(0.45)),
                corner_area(1.5, 1.0, FRAC_PI_2, conic(0.45)),
                corner_area(1.2, 0.8, FRAC_PI_2, conic(0.45)),
            ],
        ),
        (
            "two asymmetric rounds and a round",
            |f, swap, asym| {
                if swap {
                    [
                        asym(&f[1], 1.0, 2.0),
                        asym(&f[0], 1.5, 1.0),
                        Blend::Round(1.0),
                    ]
                } else {
                    [
                        asym(&f[0], 2.0, 1.0),
                        asym(&f[2], 1.0, 1.5),
                        Blend::Round(1.0),
                    ]
                }
            },
            [
                corner_area(2.0, 1.0, FRAC_PI_2, conic(0.45)),
                corner_area(1.5, 1.0, FRAC_PI_2, conic(0.45)),
                1.0 - core::f64::consts::FRAC_PI_4,
            ],
        ),
        (
            "an asymmetric round, a round and a bevel",
            |f, swap, asym| {
                let x = if swap {
                    asym(&f[1], 1.0, 2.0)
                } else {
                    asym(&f[0], 2.0, 1.0)
                };
                [
                    x,
                    Blend::Round(1.2),
                    Blend::Chamfer(ogeom_fillet::Chamfer::Symmetric(1.0)),
                ]
            },
            [
                corner_area(2.0, 1.0, FRAC_PI_2, conic(0.45)),
                (1.0 - core::f64::consts::FRAC_PI_4) * 1.44,
                0.5,
            ],
        ),
    ];
    let _ = w;
    for (name, ask, areas) in cases {
        let mut volumes = Vec::new();
        for (swap, reversed) in [(false, false), (true, false), (false, true)] {
            let mut m = Model::new();
            let (s, edges, faces) = box_corner(&mut m);
            let blends = ask(&faces, swap, &asym);
            let mut request: Vec<(Shape, Blend)> = edges.iter().cloned().zip(blends).collect();
            if reversed {
                request.reverse();
            }
            let built =
                ogeom_fillet::blend_edges(&mut m, &s, &request, ChamferCorner::Point, T).unwrap();
            let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{name}: {:?}", diagnosis.problems);
            let tangency = patch_tangency(&m, &built.shape);
            assert!(
                tangency <= 2e-4,
                "{name}: the patch meets a band at {tangency} rad"
            );
            let v = volume(&m, &built.shape);
            let lengths = [20.0, 20.0, 10.0];
            // The patch also takes the corner of the cell the bands leave,
            // no more than the cube of the widest band's reach.
            let whole: f64 = areas.iter().zip(lengths).map(|(a, l)| a * l).sum::<f64>() + 8.0;
            let widest = areas
                .iter()
                .zip(lengths)
                .map(|(a, l)| a * (l - 4.0))
                .fold(0.0, f64::max);
            assert!(
                v > 4000.0 - whole && v < 4000.0 - widest,
                "{name}: volume {v} out of ({}, {})",
                4000.0 - whole,
                4000.0 - widest
            );
            volumes.push(v.to_bits());
        }
        assert!(
            volumes.iter().all(|v| *v == volumes[0]),
            "{name}: naming the other face or reordering changes the solid"
        );
    }
}

/// Where an asymmetric round meets other rounds at a corner that is not a
/// convex corner of planes (the L-bracket's re-entrant edge, 1.5 on the
/// floor and 1 up the wall, where it meets the convex edges of its plate and
/// its upright, each round at 1), the blend engine closes the corner on one
/// patch, as it closes such a corner of rounds of one radius. Built with no
/// boolean, valid, every edge generating its strip, the same solid in either
/// order; the bracket with the fill along its 20 and less each round along
/// its 12, within the corner's own box (the setbacks' 1.5 each way).
#[test]
fn an_asymmetric_round_at_a_corner_that_is_not_convex_closes_on_one_patch() {
    let mut volumes = Vec::new();
    for flip in [false, true] {
        let mut m = Model::new();
        let (s, edge, floor, _) = l_bracket(&mut m);
        let plate = edge_where(&m, &s, |p| {
            p.y.abs() < 1e-9 && (p.z - 8.0).abs() < 1e-9 && p.x >= 8.0 - 1e-9
        });
        let upright = edge_where(&m, &s, |p| {
            p.y.abs() < 1e-9 && (p.x - 8.0).abs() < 1e-9 && p.z >= 8.0 - 1e-9
        });
        let mut request = vec![
            (
                edge,
                Blend::RoundDistances {
                    face: floor,
                    on_face: 1.5,
                    on_other: 1.0,
                    profile: RoundProfile::Elliptical,
                },
            ),
            (plate, Blend::Round(1.0)),
            (upright, Blend::Round(1.0)),
        ];
        if flip {
            request.reverse();
        }
        let (built, work) = ogeom_bool::work::counted(|| {
            ogeom_fillet::blend_edges(&mut m, &s, &request, ChamferCorner::Point, T).unwrap()
        });
        assert_eq!(work.face_pairs, 0, "built with no boolean");
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        for (e, _) in &request {
            assert!(!built.history.generated(e).is_empty(), "each strip");
        }
        volumes.push(volume(&m, &built.shape));
    }
    let spandrel = 1.0 - core::f64::consts::FRAC_PI_4;
    let want = 5120.0 + spandrel * 1.5 * 20.0 - spandrel * 24.0;
    for v in &volumes {
        assert!((v - want).abs() < 1.5 * 1.5 * 1.5, "{v} against {want}");
    }
    assert!((volumes[0] - volumes[1]).abs() < 1e-9, "either order");
}

/// The largest pair: on a box edge the top's contact runs off the top's
/// far side at a scale of the top's depth over its distance, the front's
/// at the front's height over its; on a drum's rim the wall's contact
/// runs off the wall's foot, or the cap's reaches the axis. The answer
/// builds, and a scale a hundredth past it is refused.
#[test]
fn the_largest_pair_is_where_the_first_contact_leaves_its_face() {
    let profile = RoundProfile::Conic(0.45);
    let mut m = Model::new();
    let (s, edge, top, _) = block(&mut m);
    let (on_top, on_front) = (4.0, 1.0);
    let largest =
        ogeom_fillet::largest_fillet_distances(&m, &s, &edge, &top, on_top, on_front, 100.0, T)
            .unwrap()
            .expect("a pair fits");
    assert!(
        (largest - 20.0 / on_top).abs() <= 1e-6,
        "the box edge's largest scale {largest}"
    );
    let built = ogeom_fillet::fillet_edge_distances(
        &mut m,
        &s,
        &edge,
        &top,
        largest * on_top,
        largest * on_front,
        profile,
        T,
    );
    assert!(built.is_ok(), "{:?}", built.err());
    let past = ogeom_fillet::fillet_edge_distances(
        &mut m,
        &s,
        &edge,
        &top,
        largest * on_top * 1.01,
        largest * on_front * 1.01,
        profile,
        T,
    );
    assert!(past.is_err(), "a scale past the largest is not refused");

    let mut m = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, 5.0, 4.0, T)
        .unwrap()
        .shape;
    let rim = ring_edge(&m, &drum, (0.0, 0.0, 4.0), 5.0);
    let cap = face_at(&m, &drum, Point::new(0.0, 0.0, 4.0));
    for ((on_cap, on_wall), want) in [((1.0, 2.0), 2.0), ((2.0, 0.5), 2.5)] {
        let largest = ogeom_fillet::largest_fillet_distances(
            &m, &drum, &rim, &cap, on_cap, on_wall, 100.0, T,
        )
        .unwrap()
        .expect("a pair fits");
        assert!(
            (largest - want).abs() <= 1e-6,
            "the rim's largest scale {largest}, not {want}"
        );
        let built = ogeom_fillet::fillet_edge_distances(
            &mut m,
            &drum,
            &rim,
            &cap,
            largest * on_cap,
            largest * on_wall,
            profile,
            T,
        );
        assert!(built.is_ok(), "{:?}", built.err());
    }
}

/// An asymmetric round beside a round or a bevel it does not meet: the
/// volume is the box less each section, a square less a quarter-disc or the
/// bevel's triangle and the conic corner, over the box's full width. The
/// same solid in either order, and every edge generates its band.
#[test]
fn an_asymmetric_round_beside_a_disjoint_round_or_bevel_keeps_its_history() {
    let profile = RoundProfile::Conic(0.4);
    for (other, area) in [
        (Blend::Round(2.0), 4.0 - core::f64::consts::PI),
        (Blend::Chamfer(ogeom_fillet::Chamfer::Symmetric(2.0)), 2.0),
    ] {
        let mut volumes = Vec::new();
        for reverse in [false, true] {
            let mut m = Model::new();
            let (s, edge, top, _) = block(&mut m);
            let back = edge_where(&m, &s, |p| (p.y - 20.0).abs() < 1e-9 && p.z.abs() < 1e-9);
            let asked = Blend::RoundDistances {
                face: top,
                on_face: 3.0,
                on_other: 1.5,
                profile,
            };
            let mut request = vec![(back, other.clone()), (edge, asked)];
            if reverse {
                request.reverse();
            }
            let both =
                ogeom_fillet::blend_edges(&mut m, &s, &request, ChamferCorner::Point, T).unwrap();
            let exact = 4000.0 - 20.0 * area - 20.0 * corner_area(3.0, 1.5, FRAC_PI_2, conic(0.4));
            let diagnosis = ogeom_algo::check(&m, &both.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
            let v = volume(&m, &both.shape);
            assert!((v - exact).abs() <= 1e-8, "volume {v} against {exact}");
            for (edge, _) in &request {
                assert!(both.history.is_deleted(edge));
                assert!(!both.history.generated(edge).is_empty());
            }
            volumes.push(v.to_bits());
        }
        assert_eq!(
            volumes[0], volumes[1],
            "the request's order carries no meaning"
        );
    }
}

/// A box corner closed by an asymmetric round, a round and a bevel: each of
/// the three edges is deleted and generates its band. A round on an edge
/// away from the corner takes its own section over its whole length and
/// leaves the corner as it was.
#[test]
fn a_box_corner_of_mixed_blends_keeps_each_edge_s_history() {
    let mut m = Model::new();
    let (s, edges, faces) = box_corner(&mut m);
    let back = edge_where(&m, &s, |p| (p.y - 20.0).abs() < 1e-9 && p.z.abs() < 1e-9);
    let mut request = vec![
        (
            edges[0].clone(),
            Blend::RoundDistances {
                face: faces[0].clone(),
                on_face: 2.0,
                on_other: 1.0,
                profile: RoundProfile::Conic(0.45),
            },
        ),
        (edges[1].clone(), Blend::Round(1.2)),
        (
            edges[2].clone(),
            Blend::Chamfer(ogeom_fillet::Chamfer::Symmetric(1.0)),
        ),
    ];
    let built = ogeom_fillet::blend_edges(&mut m, &s, &request, ChamferCorner::Point, T).unwrap();
    for edge in &edges {
        assert!(built.history.is_deleted(edge));
        assert!(!built.history.generated(edge).is_empty());
    }
    let corner_only = volume(&m, &built.shape);
    request.push((back.clone(), Blend::Round(2.0)));
    let with_back =
        ogeom_fillet::blend_edges(&mut m, &s, &request, ChamferCorner::Point, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &with_back.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    for edge in edges.iter().chain([&back]) {
        assert!(with_back.history.is_deleted(edge));
        assert!(!with_back.history.generated(edge).is_empty());
    }
    let exact = corner_only - 20.0 * (4.0 - core::f64::consts::PI);
    let v = volume(&m, &with_back.shape);
    assert!((v - exact).abs() <= 1e-6, "volume {v} against {exact}");
}

/// A mixed request whose asymmetric distance is not a positive finite
/// length is refused by name, whichever distance it is.
#[test]
fn a_mixed_request_with_a_distance_of_nothing_is_refused_by_name() {
    let mut m = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, 5.0, 4.0, T)
        .unwrap()
        .shape;
    let top = ring_edge(&m, &drum, (0.0, 0.0, 4.0), 5.0);
    let bottom = ring_edge(&m, &drum, (0.0, 0.0, 0.0), 5.0);
    let cap = face_at(&m, &drum, Point::new(0.0, 0.0, 4.0));
    for distance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        for on_cap in [false, true] {
            let why = ogeom_fillet::blend_edges(
                &mut m,
                &drum,
                &[
                    (bottom.clone(), Blend::Round(1.0)),
                    (
                        top.clone(),
                        Blend::RoundDistances {
                            face: cap.clone(),
                            on_face: if on_cap { distance } else { 1.5 },
                            on_other: if on_cap { 0.8 } else { distance },
                            profile: RoundProfile::Conic(0.45),
                        },
                    ),
                ],
                ChamferCorner::Point,
                T,
            )
            .unwrap_err()
            .to_string();
            assert!(
                why.contains("round set back") && why.contains("rounds nothing"),
                "{why}"
            );
        }
    }
}

/// Half a drum of radius 5 and height 4, its flat side the plane x = 0
/// through the axis: the straight edge where the flat meets the drum at
/// (0, 5, z), along the drum's ruling, the flat and the drum.
fn half_drum() -> (Model, Shape, Shape, Shape, Shape) {
    let mut m = Model::new();
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
    let wall = face_at(&m, &d, Point::new(5.0, 0.0, 2.0));
    (m, d, side, flat, wall)
}

/// The area an asymmetric round takes off, or adds to, the section of an
/// edge where a flat meets a drum of radius `r` along its ruling, the flat
/// square to the drum there: the flat's contact `d_flat` from the crease,
/// the drum's `d_drum` round its arc, a turn `beta = d_drum / r`. The
/// drum's tangent at its contact meets the flat's line `r / cos(beta)`
/// from the axis, `r (1 / cos(beta) - 1)` past the crease, at the angle
/// `pi / 2 - beta` where the crease is convex (the shoulder beyond the
/// crease, away from the flat's contact) and `pi / 2 + beta` where it is
/// concave (the shoulder between them). The conic's corner between those
/// tangents, less (convex) or with (concave) the sliver between them and
/// the drum's arc: the triangle from the axis to the shoulder and the
/// contact less the sector, `r^2 (tan(beta) - beta) / 2`.
fn ruling_section(r: f64, d_flat: f64, d_drum: f64, profile: RoundProfile, convex: bool) -> f64 {
    let beta = d_drum / r;
    let past = r / beta.cos() - r;
    let sliver = 0.5 * r * r * (beta.tan() - beta);
    let (leg, alpha) = if convex {
        (d_flat + past, FRAC_PI_2 - beta)
    } else {
        (d_flat - past, FRAC_PI_2 + beta)
    };
    let w = match profile {
        RoundProfile::Elliptical => elliptical(alpha),
        RoundProfile::Conic(rho) => conic(rho),
    };
    let corner = corner_area(leg, r * beta.tan(), alpha, w);
    if convex {
        corner - sliver
    } else {
        corner + sliver
    }
}

/// Builds `edge` of `s` naming `first` with `a` and then `second` with
/// `b`, the distances swapped; checks each valid, of `exact` volume within
/// 1e-8, its band tangent to both faces along its two straight contacts,
/// the edge consumed, and both the same solid to the bit.
fn ruled_both_ways(
    m: &mut Model,
    s: &Shape,
    edge: &Shape,
    (first, a): (&Shape, f64),
    (second, b): (&Shape, f64),
    profile: RoundProfile,
    exact: f64,
) {
    let one = ogeom_fillet::fillet_edge_distances(m, s, edge, first, a, b, profile, T).unwrap();
    let two = ogeom_fillet::fillet_edge_distances(m, s, edge, second, b, a, profile, T).unwrap();
    for built in [&one, &two] {
        let diagnosis = ogeom_algo::check(m, &built.shape, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "{profile:?}: {:?}",
            diagnosis.problems
        );
        let v = volume(m, &built.shape);
        assert!(
            (v - exact).abs() <= 1e-8,
            "{profile:?}: volume {v} against the exact {exact}"
        );
        assert!(built.history.is_deleted(edge));
        let found = contacts(m, &built.shape);
        assert_eq!(found.len(), 2, "the band meets two faces along lines");
        for (_, tangency) in &found {
            assert!(
                tangency.abs() <= 1e-7,
                "{profile:?}: the band meets a face at {tangency} rad"
            );
        }
    }
    assert_eq!(
        volume(m, &one.shape).to_bits(),
        volume(m, &two.shape).to_bits(),
        "naming the other face first"
    );
    assert_eq!(vertex_bits(m, &one.shape), vertex_bits(m, &two.shape));
}

/// The straight edge where a half drum's flat side meets its cylinder
/// runs along the drum's ruling, so every section square to it is the
/// same: the band is the shoulder conic extruded the drum's height, and
/// takes exactly its section's area over that height.
#[test]
fn an_edge_along_a_drum_s_ruling_takes_its_shoulder_conic_over_its_length() {
    let (d_flat, d_drum) = (1.0, 0.5);
    for profile in PROFILES {
        let (mut m, s, side, flat, wall) = half_drum();
        let exact =
            50.0 * core::f64::consts::PI - 4.0 * ruling_section(5.0, d_flat, d_drum, profile, true);
        ruled_both_ways(
            &mut m,
            &s,
            &side,
            (&flat, d_flat),
            (&wall, d_drum),
            profile,
            exact,
        );
    }
}

/// A drum of radius 5 standing in the corner of a slab that reaches its
/// axis, the slab's side y = 0 meeting the drum along its ruling at
/// (5, 0, z) in a concave crease: the band fills the shoulder conic's
/// section along the drum's height.
#[test]
fn a_concave_edge_along_a_drum_s_ruling_fills_its_shoulder_conic() {
    let (d_flat, d_drum) = (1.0, 0.5);
    let pi = core::f64::consts::PI;
    for profile in PROFILES {
        let mut m = Model::new();
        let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, 5.0, 4.0, T)
            .unwrap()
            .shape;
        let at = Frame::new(Point::new(0.0, -10.0, 0.0), Direction::Z, Direction::X, T).unwrap();
        let slab = ogeom_algo::make_box(&mut m, at, (10.0, 10.0, 4.0), T)
            .unwrap()
            .shape;
        let s = ogeom_bool::fuse(&mut m, &drum, &slab, T).unwrap().shape;
        let before = 75.0 * pi + 400.0;
        assert!((volume(&m, &s) - before).abs() <= 1e-8);
        let crease = edge_where(&m, &s, |p| (p.x - 5.0).abs() < 1e-9 && p.y.abs() < 1e-9);
        let side = face_at(&m, &s, Point::new(8.0, 0.0, 2.0));
        let wall = face_at(&m, &s, Point::new(0.0, 5.0, 2.0));
        let exact = before + 4.0 * ruling_section(5.0, d_flat, d_drum, profile, false);
        ruled_both_ways(
            &mut m,
            &s,
            &crease,
            (&side, d_flat),
            (&wall, d_drum),
            profile,
            exact,
        );
    }
}

/// The largest pair on seats whose legs walk along their faces: along a
/// half drum's ruling the flat's contact runs off the flat's far side, 10
/// across, at a scale of 10 over its distance; on the slab's spline
/// crease the swept wall's contact runs straight back across the wall, 2
/// deep, at 2 over its distance. The answer builds valid, and a scale a
/// hundredth past it is refused.
#[test]
fn the_largest_pair_on_a_ruling_or_a_spline_crease_is_where_a_leg_leaves_its_face() {
    let profile = RoundProfile::Elliptical;
    let (mut m, s, side, flat, _) = half_drum();
    let mut slab_m = Model::new();
    let (slab_s, crease, top, _) = slab(&mut slab_m);
    for (m, s, edge, face, (on_face, on_other), want) in [
        (&mut m, &s, &side, &flat, (1.0, 0.5), 10.0),
        (&mut slab_m, &slab_s, &crease, &top, (0.4, 0.7), 2.0 / 0.7),
    ] {
        let largest =
            ogeom_fillet::largest_fillet_distances(m, s, edge, face, on_face, on_other, 100.0, T)
                .unwrap()
                .expect("a pair fits");
        assert!(
            (largest - want).abs() <= 1e-4,
            "the largest scale {largest}, not {want}"
        );
        let built = ogeom_fillet::fillet_edge_distances(
            m,
            s,
            edge,
            face,
            largest * on_face,
            largest * on_other,
            profile,
            T,
        )
        .unwrap();
        let diagnosis = ogeom_algo::check(m, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        let past = ogeom_fillet::fillet_edge_distances(
            m,
            s,
            edge,
            face,
            largest * on_face * 1.01,
            largest * on_other * 1.01,
            profile,
            T,
        );
        assert!(past.is_err(), "a scale past the largest is not refused");
    }
}
/// How far the elliptical section of the half drum's convex ruling stands
/// out of the drum at its worst, asked 0.1 along the flat and 1 round the
/// drum per unit of `scale`: the conic from the flat's contact to the
/// drum's, about the shoulder `5 / cos(beta)` up the flat's line, of
/// weight `sin(alpha / 2)` at the shoulder's angle `alpha = pi / 2 - beta`,
/// sampled densely and measured from the axis against the radius. Its two
/// contacts lie on the drum's circle or inside it, so a positive answer is
/// the arc crossing the drum between them.
fn ruling_section_outside(scale: f64) -> f64 {
    let r = 5.0;
    let beta = scale / r;
    let p0 = (0.0, 0.1f64.mul_add(-scale, r));
    let c = (0.0, r / beta.cos());
    let p2 = (r * beta.sin(), r * beta.cos());
    let w = elliptical(FRAC_PI_2 - beta);
    let mut worst = f64::NEG_INFINITY;
    for i in 1..20_000 {
        let t = f64::from(i) / 20_000.0;
        let (a, b, e) = ((1.0 - t) * (1.0 - t), 2.0 * t * (1.0 - t) * w, t * t);
        let sum = a + b + e;
        let x = (a * p0.0 + b * c.0 + e * p2.0) / sum;
        let y = (a * p0.1 + b * c.1 + e * p2.1) / sum;
        worst = worst.max(x.hypot(y) - r);
    }
    worst
}

/// The scale at which the section of [`ruling_section_outside`] first
/// reaches the drum between its contacts, bisected.
fn ruling_section_meets_drum() -> f64 {
    let (mut lo, mut hi) = (1.0, 1.8);
    assert!(ruling_section_outside(lo) <= 0.0 && ruling_section_outside(hi) > 0.0);
    for _ in 0..60 {
        let mid = f64::midpoint(lo, hi);
        if ruling_section_outside(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    lo
}

/// The half drum's convex ruling asked 0.1 along the flat and 1 round
/// the drum per unit of scale, the drum's leg turning toward a quarter
/// turn. The conic standing on the shoulder bulges toward the drum as the
/// turn grows, and from the scale [`ruling_section_meets_drum`] (about
/// 1.6914) it crosses the drum between the contacts: no band of that
/// section is tangent to the drum there without leaving it. Up to that
/// scale the band builds valid, of the section's closed-form area over
/// the drum's height; past it the round is refused by name, never built
/// bulging out of the drum.
#[test]
fn a_convex_ruling_builds_until_its_section_reaches_the_drum() {
    let meets = ruling_section_meets_drum();
    let profile = RoundProfile::Elliptical;
    for scale in [1.0, 1.5, 1.68, 1.69, meets * (1.0 - 2e-4)] {
        let (mut m, s, side, flat, _) = half_drum();
        let built = ogeom_fillet::fillet_edge_distances(
            &mut m,
            &s,
            &side,
            &flat,
            0.1 * scale,
            scale,
            profile,
            T,
        )
        .unwrap_or_else(|e| panic!("scale {scale}: {e}"));
        let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "scale {scale}: {:?}",
            diagnosis.problems
        );
        let exact = 50.0 * core::f64::consts::PI
            - 4.0 * ruling_section(5.0, 0.1 * scale, scale, profile, true);
        let v = volume(&m, &built.shape);
        assert!(
            (v - exact).abs() <= 1e-8,
            "scale {scale}: volume {v} against the exact {exact}"
        );
        for (_, tangency) in &contacts(&m, &built.shape) {
            assert!(
                tangency.abs() <= 1e-7,
                "scale {scale}: meets at {tangency} rad"
            );
        }
    }
    for scale in [meets * 1.001, 1.8, 1.9085, 3.0, 5.0, 7.8] {
        let (mut m, s, side, flat, _) = half_drum();
        let message = ogeom_fillet::fillet_edge_distances(
            &mut m,
            &s,
            &side,
            &flat,
            0.1 * scale,
            scale,
            profile,
            T,
        )
        .expect_err("refused")
        .to_string();
        assert!(
            message.contains("crosses the drum"),
            "scale {scale}: {message}"
        );
    }
}

/// The largest pair on the half drum's convex ruling asked 0.1 along the
/// flat and 1 round the drum: the scale at which the section reaches the
/// drum ([`ruling_section_meets_drum`]), not the shoulder's quarter turn
/// nor where a boolean gives up. The answer builds valid, and a scale a
/// hundredth past it is refused.
#[test]
fn the_largest_pair_on_a_convex_ruling_stops_where_its_section_reaches_the_drum() {
    let want = ruling_section_meets_drum();
    let (mut m, s, side, flat, _) = half_drum();
    let largest = ogeom_fillet::largest_fillet_distances(&m, &s, &side, &flat, 0.1, 1.0, 100.0, T)
        .unwrap()
        .expect("a pair fits");
    assert!(
        (largest - want).abs() <= 1e-4 * want,
        "the largest scale {largest}, not {want}"
    );
    let profile = RoundProfile::Elliptical;
    let built = ogeom_fillet::fillet_edge_distances(
        &mut m,
        &s,
        &side,
        &flat,
        0.1 * largest,
        largest,
        profile,
        T,
    )
    .unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let past = ogeom_fillet::fillet_edge_distances(
        &mut m,
        &s,
        &side,
        &flat,
        0.1 * largest * 1.01,
        largest * 1.01,
        profile,
        T,
    );
    assert!(past.is_err(), "a scale past the largest is not refused");
}

/// How far a square corner's conic section of legs `a` along one face and
/// `b` along the other, of weight `w`, stands off the first face at `u`
/// along it from the crease: the corner's depth there, 0 past `a`. Read
/// off the rational arc from the first face's contact to the second's,
/// bisected on its parameter, along which `u` falls.
fn conic_depth(a: f64, b: f64, w: f64, u: f64) -> f64 {
    if u >= a {
        return 0.0;
    }
    if u <= 0.0 {
        return b;
    }
    let at = |t: f64| {
        let (p, q, r) = ((1.0 - t) * (1.0 - t), 2.0 * t * (1.0 - t) * w, t * t);
        let sum = p + q + r;
        (p * a / sum, r * b / sum)
    };
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if at(mid).0 > u {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    at(0.5 * (lo + hi)).1
}

/// Gauss-Legendre over `[0, top]` in 400 pieces, the variable laid as
/// `top s^2` so the square-root ends of a conic's depth integrate smooth.
fn integrate(top: f64, f: impl Fn(f64) -> f64) -> f64 {
    let nodes = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
    ];
    let mut sum = 0.0;
    for k in 0..400 {
        let (a, b) = (f64::from(k) / 400.0, f64::from(k + 1) / 400.0);
        for (x, wt) in nodes {
            let s = (b - a).mul_add(x, a + b) * 0.5;
            sum += wt * (b - a) * 0.5 * f(top * s * s) * 2.0 * top * s;
        }
    }
    sum
}

/// Builds `request` on `s` in its order and reversed: each valid, every
/// edge generating its band, the same solid either way; its volume.
fn either_order(m: &mut Model, s: &Shape, request: &[(Shape, Blend)]) -> f64 {
    let mut volumes = Vec::new();
    let mut reversed = request.to_vec();
    reversed.reverse();
    for asked in [request.to_vec(), reversed] {
        let built = ogeom_fillet::blend_edges(m, s, &asked, ChamferCorner::Point, T).unwrap();
        let found = ogeom_algo::check(m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        for (e, _) in &asked {
            assert!(!built.history.generated(e).is_empty(), "each band");
        }
        volumes.push(volume(m, &built.shape));
    }
    assert!(
        (volumes[0] - volumes[1]).abs() <= 1e-9 * volumes[0],
        "either order: {volumes:?}"
    );
    volumes[0]
}

/// Two of a box corner's three edges blended, the third left sharp: the
/// two bands run to their end faces and meet by intersection, so the box
/// loses each band's corner along its whole edge and gets back once what
/// both take near the corner. Beside a round of the upright edge, that is
/// a column over each station of the top's left face: the asymmetric
/// band's depth there times the round's width. Beside an asymmetric round
/// of the top's left edge, each point of the top under both takes the
/// shallower band's depth.
#[test]
fn bands_where_the_third_edge_is_left_sharp_meet_by_intersection() {
    let w = conic(0.4);
    let asym = |face: &Shape, a: f64, b: f64| Blend::RoundDistances {
        face: face.clone(),
        on_face: a,
        on_other: b,
        profile: RoundProfile::Conic(0.4),
    };
    let circle = FRAC_1_SQRT_2;
    let (a_top, b_front, r) = (3.0, 1.5, 1.0);
    let x_corner = corner_area(a_top, b_front, FRAC_PI_2, w);

    // With a round of the upright edge, the top's left edge sharp.
    let mut m = Model::new();
    let (s, [x, _, z], [top, ..]) = box_corner(&mut m);
    let both = integrate(r, |y| {
        conic_depth(a_top, b_front, w, y) * conic_depth(r, r, circle, y)
    });
    let exact =
        4000.0 - 20.0 * x_corner - 10.0 * r * r * (1.0 - core::f64::consts::FRAC_PI_4) + both;
    let v = either_order(
        &mut m,
        &s,
        &[(x, asym(&top, a_top, b_front)), (z, Blend::Round(r))],
    );
    assert!((v - exact).abs() <= 1e-7, "volume {v} against {exact}");

    // With an asymmetric round of the top's left edge, the upright sharp:
    // 1 along the top, 2 down the left face.
    let mut m = Model::new();
    let (s, [x, y, _], [top, _, left]) = box_corner(&mut m);
    let (a_left, b_top) = (2.0, 1.0);
    let both = integrate(b_top, |u| {
        let k = conic_depth(b_top, a_left, w, u);
        integrate(a_top, |t| conic_depth(a_top, b_front, w, t).min(k))
    });
    let exact = 4000.0 - 20.0 * x_corner - 20.0 * corner_area(a_left, b_top, FRAC_PI_2, w) + both;
    let v = either_order(
        &mut m,
        &s,
        &[
            (x, asym(&top, a_top, b_front)),
            (y, asym(&left, a_left, b_top)),
        ],
    );
    assert!((v - exact).abs() <= 1e-5, "volume {v} against {exact}");
}

/// The half of a drum of radius 5 and height 4 past x = 0, its seam
/// behind the axis so its top rim is one arc: the half drum, the top's arc
/// and straight edges, the top, the flat and the drum's face.
fn d_shape(m: &mut Model) -> (Shape, Shape, Shape, [Shape; 3]) {
    let back = Frame::new(
        Point::new(0.0, 0.0, 0.0),
        Direction::Z,
        Direction::new(-X, T).unwrap(),
        T,
    )
    .unwrap();
    let drum = ogeom_algo::make_cylinder(m, back, 5.0, 4.0, T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(-6.0, -6.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let half = ogeom_algo::make_box(m, at, (6.0, 12.0, 6.0), T)
        .unwrap()
        .shape;
    let d = ogeom_bool::cut(m, &drum, &half, T).unwrap().shape;
    let mid = |m: &Model, p: Point| {
        use ogeom_geom::Curve3d as _;
        explore_unique(m, &d, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                ogeom_algo::placed_edge_curve(m, e, T)
                    .unwrap()
                    .is_some_and(|(c, (lo, hi))| {
                        c.point_at(f64::midpoint(lo, hi), T).unwrap().distance(p) < 1e-6
                    })
            })
            .expect("the edge")
    };
    let arc = mid(m, Point::new(5.0, 0.0, 4.0));
    let straight = mid(m, Point::new(0.0, 0.0, 4.0));
    let top = face_at(m, &d, Point::new(2.0, 0.0, 4.0));
    let flat = face_at(m, &d, Point::new(0.0, 0.0, 2.0));
    let wall = face_at(m, &d, Point::new(5.0, 0.0, 2.0));
    (d, arc, straight, [top, flat, wall])
}

/// A curved edge at a vertex left sharp: on the D-shaped half drum, an
/// asymmetric round of the top's arc and a round or asymmetric round of
/// its straight edge, the uprights at the arc's ends left sharp. The arc's
/// band is its conic revolved over the arc's half turn, ending in the flat
/// through the axis; the straight band runs to the drum, square to it
/// there. Both hang from the top, so the solid loses at each point of the
/// top the deeper of the two bands' depths there: the straight band's at
/// its distance from the flat, the arc's at its distance in from the rim.
/// That is integrated over the half disc, split where either band's
/// contact or their crossing runs. The same with a circular round on the
/// arc, its band the torus. The same solid in either order and naming the
/// other faces. A bevel on the straight edge, its prism ending at the
/// drum, is refused by name.
#[test]
fn an_asymmetric_round_of_a_rim_s_arc_beside_an_edge_left_sharp_meets_the_other_band_by_intersection()
 {
    let w = conic(0.4);
    let asym = |face: &Shape, a: f64, b: f64| Blend::RoundDistances {
        face: face.clone(),
        on_face: a,
        on_other: b,
        profile: RoundProfile::Conic(0.4),
    };
    let (cap, rim_wall) = (1.0, 0.6);
    let arc_depth = |rho: f64| conic_depth(cap, rim_wall, w, 5.0 - rho);
    let taken_by = |straight: &dyn Fn(f64) -> f64,
                    reach: f64,
                    arc_depth: &dyn Fn(f64) -> f64,
                    cap: f64| {
        // Over y >= 0 at each x, the arc's depth grows toward the rim.
        let across = |x: f64| {
            let top = x.mul_add(-x, 25.0).max(0.0).sqrt();
            let h = straight(x);
            let mut cuts = vec![0.0, top];
            let inner = 5.0 - cap;
            if x < inner {
                cuts.push(inner.mul_add(inner, -x * x).sqrt());
            }
            let gap = |y: f64| arc_depth(x.hypot(y).min(5.0)) - h;
            if h > 0.0 && gap(0.0) < 0.0 && gap(top) > 0.0 {
                let (mut lo, mut hi) = (0.0, top);
                for _ in 0..80 {
                    let mid = f64::midpoint(lo, hi);
                    if gap(mid) > 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                    }
                }
                cuts.push(f64::midpoint(lo, hi));
            }
            cuts.sort_by(f64::total_cmp);
            cuts.windows(2)
                .map(|p| integrate_between(p[0], p[1], |y| h.max(arc_depth(x.hypot(y).min(5.0)))))
                .sum::<f64>()
        };
        let mut xs = [0.0, reach.min(5.0), 5.0 - cap, 5.0];
        xs.sort_by(f64::total_cmp);
        2.0 * xs
            .windows(2)
            .map(|p| integrate_between(p[0], p[1], across))
            .sum::<f64>()
    };
    let taken =
        |straight: &dyn Fn(f64) -> f64, reach: f64| taken_by(straight, reach, &arc_depth, cap);
    let solid = 50.0 * core::f64::consts::PI;
    let circle = FRAC_1_SQRT_2;

    // A round of radius 0.8 on the straight edge.
    let mut m = Model::new();
    let (d, arc, straight, [top, flat, wall]) = d_shape(&mut m);
    let exact = solid - taken(&|x| conic_depth(0.8, 0.8, circle, x), 0.8);
    let v = either_order(
        &mut m,
        &d,
        &[
            (arc.clone(), asym(&top, cap, rim_wall)),
            (straight.clone(), Blend::Round(0.8)),
        ],
    );
    assert!(
        (v - exact).abs() <= 1e-7,
        "round: volume {v} against {exact}"
    );
    let again = ogeom_fillet::blend_edges(
        &mut m,
        &d,
        &[
            (arc.clone(), asym(&wall, rim_wall, cap)),
            (straight.clone(), Blend::Round(0.8)),
        ],
        ChamferCorner::Point,
        T,
    )
    .unwrap();
    assert_eq!(
        volume(&m, &again.shape).to_bits(),
        v.to_bits(),
        "naming the wall"
    );

    // An asymmetric round of 0.9 along the top and 0.5 down the flat.
    let exact = solid - taken(&|x| conic_depth(0.9, 0.5, w, x), 0.9);
    let v = either_order(
        &mut m,
        &d,
        &[
            (arc.clone(), asym(&top, cap, rim_wall)),
            (straight.clone(), asym(&flat, 0.5, 0.9)),
        ],
    );
    assert!(
        (v - exact).abs() <= 1e-7,
        "asymmetric: volume {v} against {exact}"
    );

    let refused = |m: &mut Model, request: &[(Shape, Blend)], words: &str| {
        let message = ogeom_fillet::blend_edges(m, &d, request, ChamferCorner::Point, T)
            .expect_err("refused")
            .to_string();
        assert!(message.contains(words), "{message}");
    };
    // A circular round of 0.8 on the arc, its band the torus.
    let round_arc = |rho: f64| conic_depth(0.8, 0.8, circle, 5.0 - rho);
    let exact = solid - taken_by(&|x| conic_depth(0.9, 0.5, w, x), 0.9, &round_arc, 0.8);
    let v = either_order(
        &mut m,
        &d,
        &[
            (arc.clone(), Blend::Round(0.8)),
            (straight.clone(), asym(&top, 0.9, 0.5)),
        ],
    );
    assert!(
        (v - exact).abs() <= 1e-7,
        "round arc: volume {v} against {exact}"
    );
    refused(
        &mut m,
        &[
            (arc, asym(&top, cap, rim_wall)),
            (straight, bevel(&flat, 0.5, 0.9)),
        ],
        "its prism stops square there",
    );
}

/// A 20 x 20 x 10 block with the quarter above z = 5 past x = 10 and
/// y = 10 taken out: the pocket's inside corner at (10, 10, 5), its three
/// concave edges along x, y and z, the floor and the walls x = 10 and
/// y = 10.
fn pocket(m: &mut Model) -> (Shape, [Shape; 3], [Shape; 3]) {
    let block = ogeom_algo::make_box(m, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(10.0, 10.0, 5.0), Direction::Z, Direction::X, T).unwrap();
    let notch = ogeom_algo::make_box(m, at, (20.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(m, &block, &notch, T).unwrap().shape;
    let ex = edge_where(m, &s, |p| {
        (p.y - 10.0).abs() < 1e-9 && (p.z - 5.0).abs() < 1e-9 && p.x >= 10.0 - 1e-9
    });
    let ey = edge_where(m, &s, |p| {
        (p.x - 10.0).abs() < 1e-9 && (p.z - 5.0).abs() < 1e-9 && p.y >= 10.0 - 1e-9
    });
    let ez = edge_where(m, &s, |p| {
        (p.x - 10.0).abs() < 1e-9 && (p.y - 10.0).abs() < 1e-9 && p.z >= 5.0 - 1e-9
    });
    let floor = face_at(m, &s, Point::new(15.0, 15.0, 5.0));
    let wall_x = face_at(m, &s, Point::new(10.0, 15.0, 7.5));
    let wall_y = face_at(m, &s, Point::new(15.0, 10.0, 7.5));
    (s, [ex, ey, ez], [floor, wall_x, wall_y])
}

/// Two of a pocket's inside corner's three concave edges filled, the
/// upright left sharp: each fill runs to its end faces and the two meet by
/// intersection, so the block gains each fill's corner along its edge less
/// once what both fill near the corner, the shallower fill's depth over
/// each point of the floor under both.
#[test]
fn fills_where_the_third_edge_is_left_sharp_meet_by_intersection() {
    let w = conic(0.4);
    let mut m = Model::new();
    let (s, [ex, ey, _], [floor, ..]) = pocket(&mut m);
    let (a_floor, b_wall, r) = (1.5, 1.0, 1.0);
    let request = [
        (
            ex,
            Blend::RoundDistances {
                face: floor,
                on_face: a_floor,
                on_other: b_wall,
                profile: RoundProfile::Conic(0.4),
            },
        ),
        (ey, Blend::Round(r)),
    ];
    let both = integrate(r, |u| {
        let k = conic_depth(r, r, FRAC_1_SQRT_2, u);
        integrate(a_floor, |t| conic_depth(a_floor, b_wall, w, t).min(k))
    });
    let exact = 3500.0
        + 10.0 * corner_area(a_floor, b_wall, FRAC_PI_2, w)
        + 10.0 * r * r * (1.0 - core::f64::consts::FRAC_PI_4)
        - both;
    let v = either_order(&mut m, &s, &request);
    assert!((v - exact).abs() <= 1e-5, "volume {v} against {exact}");
}

/// Every edge at a pocket's inside corner blended, an asymmetric round
/// among them, or at a mixed corner where the concave upright meets two
/// convex top edges: each band stops at its setback and one patch closes
/// the corner, tangent to every band it meets. Past the largest distance
/// `m` from the corner the bands are apart and each takes or adds its own
/// section along its edge; the volume less those is what the corner's cube
/// of side `m` holds, no more than the cube and, where every band fills,
/// no less than the largest band's share of it. The same solid in either
/// order, and naming the other face.
#[test]
fn asymmetric_rounds_close_concave_and_mixed_corners_on_one_patch() {
    let w = conic(0.4);
    let asym = |face: &Shape, a: f64, b: f64| Blend::RoundDistances {
        face: face.clone(),
        on_face: a,
        on_other: b,
        profile: RoundProfile::Conic(0.4),
    };
    let spandrel = 1.0 - core::f64::consts::FRAC_PI_4;
    let ring = corner_area(1.5, 1.0, FRAC_PI_2, w);
    let side = 1.5_f64;
    let cube = side.powi(3);
    let checks = |m: &mut Model, s: &Shape, request: &[(Shape, Blend)]| -> f64 {
        let v = either_order(m, s, request);
        let built = ogeom_fillet::blend_edges(m, s, request, ChamferCorner::Point, T).unwrap();
        let tangency = patch_tangency(m, &built.shape);
        assert!(tangency <= 2e-4, "the patch meets a band at {tangency} rad");
        v
    };

    // The inside corner: an asymmetric round and two rounds, then three
    // asymmetric rounds; each fills.
    let mut m = Model::new();
    let (s, [ex, ey, ez], [floor, wall_x, wall_y]) = pocket(&mut m);
    let v = checks(
        &mut m,
        &s,
        &[
            (ex.clone(), asym(&floor, 1.5, 1.0)),
            (ey.clone(), Blend::Round(1.0)),
            (ez.clone(), Blend::Round(1.0)),
        ],
    );
    let near = v - 3500.0 - (ring * 8.5 + spandrel * 8.5 + spandrel * 3.5);
    assert!(
        ring * side <= near && near <= cube,
        "{near} at the corner against [{}, {cube}]",
        ring * side
    );
    let ez_corner = corner_area(1.2, 0.8, FRAC_PI_2, w);
    let three = [
        (ex.clone(), asym(&floor, 1.5, 1.0)),
        (ey.clone(), asym(&floor, 1.5, 1.0)),
        (ez.clone(), asym(&wall_x, 1.2, 0.8)),
    ];
    let v = checks(&mut m, &s, &three);
    let near = v - 3500.0 - (2.0 * ring * 8.5 + ez_corner * 3.5);
    assert!(
        ring * side <= near && near <= cube,
        "{near} at the corner against [{}, {cube}]",
        ring * side
    );
    // Naming the upright's other wall, the distances swapped.
    let swapped = [
        three[0].clone(),
        three[1].clone(),
        (ez.clone(), asym(&wall_y, 0.8, 1.2)),
    ];
    let again = ogeom_fillet::blend_edges(&mut m, &s, &swapped, ChamferCorner::Point, T).unwrap();
    assert_eq!(
        volume(&m, &again.shape).to_bits(),
        v.to_bits(),
        "naming the other face"
    );

    // The mixed corner on top: the concave upright filled by an
    // asymmetric round, the two convex top edges cut by rounds.
    let tx = edge_where(&m, &s, |p| {
        (p.y - 10.0).abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9 && p.x >= 10.0 - 1e-9
    });
    let ty = edge_where(&m, &s, |p| {
        (p.x - 10.0).abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9 && p.y >= 10.0 - 1e-9
    });
    let v = checks(
        &mut m,
        &s,
        &[
            (ez, asym(&wall_x, 1.5, 1.0)),
            (tx, Blend::Round(1.0)),
            (ty, Blend::Round(1.0)),
        ],
    );
    let near = v - 3500.0 - (ring * 3.5 - 2.0 * spandrel * 8.5);
    assert!(
        near.abs() <= cube,
        "{near} at the corner, past the cube {cube}"
    );
}
/// A drum's top rim asked with distances following closed laws through
/// two points each, 1 on the cap and 0.5 on the wall at the rim's start,
/// 1.5 and 0.8 half way round. With two points a closed law is level at
/// both, each half between them the cubic `a + (b - a)(3 t^2 - 2 t^3)`.
/// Every section of the band stands in a plane through the axis, so the
/// band takes `int A(s) (R - c(s)) 2 pi ds` over the turn, `A` the conic
/// corner of the station's distances and `c` its centroid's depth into
/// the cap, whichever station the rim starts at.
#[test]
fn distances_varying_round_a_rim_take_the_integral_of_their_corners() {
    use core::f64::consts::{PI, TAU};
    use ogeom_fillet::RadiusLaw;
    let (radius, height) = (5.0, 4.0);
    let on_cap = RadiusLaw::Points(vec![(0.0, 1.0), (0.5, 1.5)]);
    let on_wall = RadiusLaw::Points(vec![(0.0, 0.5), (0.5, 0.8)]);
    let smooth = |a: f64, b: f64, s: f64| {
        let t = if s <= 0.5 { s / 0.5 } else { (1.0 - s) / 0.5 };
        (b - a).mul_add(t * t * (3.0 - 2.0 * t), a)
    };
    for profile in PROFILES {
        let w = weight_of(profile);
        let taken = 2.0
            * integrate(0.5, |s| {
                let (d_cap, d_wall) = (smooth(1.0, 1.5, s), smooth(0.5, 0.8, s));
                let (area, centroid) = square_corner(d_cap, d_wall, w);
                area * (radius - centroid[0]) * TAU
            });
        let mut m = Model::new();
        let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, radius, height, T)
            .unwrap()
            .shape;
        let rim = ring_edge(&m, &drum, (0.0, 0.0, height), radius);
        let cap = face_at(&m, &drum, Point::new(0.0, 0.0, height));
        let wall = face_at(&m, &drum, Point::new(0.0, radius, 1.0));
        let exact = PI * radius * radius * height - taken;
        let mut volumes = Vec::new();
        for (face, a, b) in [(&cap, &on_cap, &on_wall), (&wall, &on_wall, &on_cap)] {
            let built = ogeom_fillet::fillet_edge_distances_variable(
                &mut m, &drum, &rim, face, a, b, profile, T,
            )
            .unwrap();
            let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
            assert!(
                diagnosis.is_valid(),
                "{profile:?}: {:?}",
                diagnosis.problems
            );
            assert!(built.history.is_deleted(&rim));
            let v = volume(&m, &built.shape);
            assert!(
                (v - exact).abs() <= 1e-4 * taken,
                "{profile:?}: volume {v} against {exact}"
            );
            volumes.push(v);
        }
        assert!(
            (volumes[0] - volumes[1]).abs() <= 1e-9 * volumes[0],
            "naming the other face first"
        );
    }
}

/// The slab's spline crease asked with distances rising linearly along its
/// length, 0.3 to 0.6 on the flat side and 0.5 to 0.9 on the swept wall,
/// square to each other all along. Each station's section is the square
/// conic corner of its distances, standing in the plane square to the
/// crease, and the band takes `int (A - k M) d sigma` over the crease's
/// arc length: `A` the corner's area, `M` its first moment about the
/// crease along the flat side's leg, and `k` the crease's curvature
/// toward that leg.
#[test]
fn distances_varying_along_a_spline_crease_take_the_integral_of_their_corners() {
    use ogeom_fillet::RadiusLaw;
    let flat_law = RadiusLaw::Linear {
        start: 0.3,
        end: 0.6,
    };
    let wall_law = RadiusLaw::Linear {
        start: 0.5,
        end: 0.9,
    };
    // The spline's Bezier points in (x, z), its derivatives, and its
    // arc length by Gauss-Legendre.
    let p = [(0.0, 0.0), (3.0, 0.0), (7.0, -1.0), (10.0, -1.0)];
    let d1 = |t: f64| {
        let u = 1.0 - t;
        let c = [3.0 * u * u, 6.0 * u * t, 3.0 * t * t];
        let x = c[0] * (p[1].0 - p[0].0) + c[1] * (p[2].0 - p[1].0) + c[2] * (p[3].0 - p[2].0);
        let z = c[0] * (p[1].1 - p[0].1) + c[1] * (p[2].1 - p[1].1) + c[2] * (p[3].1 - p[2].1);
        (x, z)
    };
    let d2 = |t: f64| {
        let c = [6.0 * (1.0 - t), 6.0 * t];
        let x = c[0] * (p[2].0 - 2.0 * p[1].0 + p[0].0) + c[1] * (p[3].0 - 2.0 * p[2].0 + p[1].0);
        let z = c[0] * (p[2].1 - 2.0 * p[1].1 + p[0].1) + c[1] * (p[3].1 - 2.0 * p[2].1 + p[1].1);
        (x, z)
    };
    let speed = |t: f64| {
        let (x, z) = d1(t);
        x.hypot(z)
    };
    let gauss = |a: f64, b: f64, f: &dyn Fn(f64) -> f64| {
        let nodes = [
            (-0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
            (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
            (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
            (0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
        ];
        nodes
            .iter()
            .map(|(x, wt)| wt * (b - a) * 0.5 * f((b - a).mul_add(*x, a + b) * 0.5))
            .sum::<f64>()
    };
    let pieces = 512;
    let mut start = vec![0.0];
    for k in 0..pieces {
        let (a, b) = (
            f64::from(k) / f64::from(pieces),
            f64::from(k + 1) / f64::from(pieces),
        );
        start.push(start[k as usize] + gauss(a, b, &speed));
    }
    let length = start[pieces as usize];
    for profile in PROFILES {
        let w = weight_of(profile);
        let mut m = Model::new();
        let (solid, crease, flat, wall) = slab(&mut m);
        // The law runs from the edge's start.
        let (a, _) = ogeom_algo::edge_vertices(&m, &crease).unwrap().unwrap();
        let from_left = point_of(&m, &a).x < 5.0;
        let mut taken = 0.0;
        for k in 0..pieces {
            let (a, b) = (
                f64::from(k) / f64::from(pieces),
                f64::from(k + 1) / f64::from(pieces),
            );
            taken += gauss(a, b, &|t| {
                let sigma = start[k as usize] + gauss(a, t, &speed);
                let s = if from_left {
                    sigma / length
                } else {
                    1.0 - sigma / length
                };
                let (d_flat, d_wall) = (0.3 + 0.3 * s, 0.5 + 0.4 * s);
                let (area, centroid) = square_corner(d_flat, d_wall, w);
                // The unit tangent, the flat's leg down from the crease,
                // and the crease's turn toward it per unit length.
                let v = speed(t);
                let (tx, tz) = (d1(t).0 / v, d1(t).1 / v);
                let (ax, az) = d2(t);
                let along = ax * tx + az * tz;
                let (nx, nz) = ((ax - along * tx) / (v * v), (az - along * tz) / (v * v));
                let toward = nx * tz - nz * tx;
                (area - toward * area * centroid[0]) * v
            });
        }
        let before = volume(&m, &solid);
        let exact = before - taken;
        let mut volumes = Vec::new();
        for (face, a, b) in [(&flat, &flat_law, &wall_law), (&wall, &wall_law, &flat_law)] {
            let built = ogeom_fillet::fillet_edge_distances_variable(
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
                (v - exact).abs() <= 1e-4 * taken,
                "{profile:?}: volume {v} against {exact}, taking {taken}"
            );
            volumes.push(v);
        }
        assert!(
            (volumes[0] - volumes[1]).abs() <= 1e-9 * volumes[0],
            "naming the other face first"
        );
    }
}

/// A 40 x 20 x 5 slab on the leaning block, its front y = 5 standing on
/// the block's top: the concave edge along (x, 5, 10) from x = 0 to where
/// the block's end crosses it, at 20 + lean / 4. With the top face and
/// the slab's front.
fn slab_on_leaning_block(m: &mut Model, lean: f64) -> (Shape, Shape, Shape, Shape) {
    let leaning = leaning_block(m, lean);
    let slab = ogeom_algo::make_box(
        m,
        Frame::new(Point::new(-5.0, 5.0, 10.0), Direction::Z, Direction::X, T).unwrap(),
        (40.0, 20.0, 5.0),
        T,
    )
    .unwrap()
    .shape;
    let s = ogeom_bool::fuse(m, &leaning, &slab, T).unwrap().shape;
    let edge = edge_where(m, &s, |p| {
        (p.y - 5.0).abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9
    });
    let top = face_at(m, &s, Point::new(10.0, 2.5, 10.0));
    let wall = face_at(m, &s, Point::new(10.0, 5.0, 12.5));
    (s, edge, top, wall)
}

/// The worst angle between a band and the top (z = 10) or the slab's
/// front (y = 5) along the curves it meets them in. Faces on planes and
/// cylinders are the solid's own.
fn tangent_to_top_and_wall(m: &Model, s: &Shape) -> f64 {
    let mut worst = 0.0_f64;
    for f in explore_unique(m, s, ShapeType::Face).unwrap() {
        let data = m.node(&f).unwrap().data().as_face().unwrap();
        if matches!(
            m.geometry().surface(data.surface),
            Some(ogeom_geom::SurfaceGeometry::Plane(_) | ogeom_geom::SurfaceGeometry::Cylinder(_))
        ) {
            continue;
        }
        for c in ogeom_fillet::analyse_blend(m, s, &f, 17, T).unwrap() {
            let (curve, (lo, hi)) = ogeom_algo::placed_edge_curve(m, &c.edge, T)
                .unwrap()
                .unwrap();
            use ogeom_geom::Curve3d as _;
            let p = curve.point_at(f64::midpoint(lo, hi), T).unwrap();
            let on_top = (p.z - 10.0).abs() <= 1e-7 && p.y < 5.0 - 1e-7;
            let on_wall = (p.y - 5.0).abs() <= 1e-7 && p.z > 10.0 + 1e-7;
            if on_top || on_wall {
                worst = worst.max(c.tangency_error);
            }
        }
    }
    worst
}

/// A concave edge ending at a face not square to it: the fill ends in that
/// face's plane, as a convex band's end is cut by the face it runs out
/// across. On the slab on the leaning block the end through
/// (20 + lean / 4, 5, 10) stands back from the crease by `k u`, `k =
/// lean / 20`, at a depth `u` into the top. Where it leans out (`k < 0`)
/// the fill runs on past the edge's end with its last section until the
/// plane cuts it, each point of that section `-k u` further; where it
/// leans in (`k > 0`) each section keeps only what lies within `(L - x) /
/// k` of the crease. With constant distances that is the corner `A` along
/// `L - k c`, `c` its centroid's depth into the top; with the laws 1 to 2
/// on the top and 0.5 to 1.5 on the slab, each station's corner over the
/// edge, the end's pieces integrated as above. The same solid whichever
/// face is named.
#[test]
fn a_concave_edge_ending_at_a_face_not_square_fills_up_to_its_plane() {
    use ogeom_fillet::RadiusLaw;
    let rho = 0.4;
    let w = conic(rho);
    let profile = RoundProfile::Conic(rho);
    let kept = corner_area(1.0, 1.0, FRAC_PI_2, w) / 0.5;
    for lean in [5.0, -5.0] {
        let k = lean / 20.0;
        let length = 20.0 + lean / 4.0;
        let solid = 100.0 * (40.0 + lean) + 4000.0;

        // Constant distances.
        let mut m = Model::new();
        let (s, edge, top, wall) = slab_on_leaning_block(&mut m, lean);
        assert!((volume(&m, &s) - solid).abs() <= 1e-8);
        let (area, c) = square_corner(1.5, 1.0, w);
        let exact = solid + area * (length - k * c[0]);
        let mut bits = Vec::new();
        for (face, a, b) in [(&top, 1.5, 1.0), (&wall, 1.0, 1.5)] {
            let built =
                ogeom_fillet::fillet_edge_distances(&mut m, &s, &edge, face, a, b, profile, T)
                    .unwrap();
            let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{lean}: {:?}", diagnosis.problems);
            let v = volume(&m, &built.shape);
            assert!(
                (v - exact).abs() <= 1e-8,
                "{lean}: constant, volume {v} against {exact}"
            );
            let tangency = tangent_to_top_and_wall(&m, &built.shape);
            assert!(
                tangency <= 1e-7,
                "{lean}: the fill meets a face at {tangency}"
            );
            bits.push(v.to_bits());
        }
        assert_eq!(bits[0], bits[1], "{lean}: naming the other face first");

        // Distances varying.
        let on_top = RadiusLaw::Linear {
            start: 1.0,
            end: 2.0,
        };
        let on_wall = RadiusLaw::Linear {
            start: 0.5,
            end: 1.5,
        };
        // L int_0^t (1 + s)(0.5 + s) ds.
        let product = |t: f64| length * t * (0.5 + t * (0.75 + t / 3.0));
        let added = if k < 0.0 {
            let (area, c) = square_corner(2.0, 1.5, w);
            0.5 * kept * product(1.0) - k * area * c[0]
        } else {
            // Whole sections up to where the plane reaches the top's
            // contact, then what of each lies within (L - x) / k.
            let whole = (length - k) / (1.0 + k / length);
            0.5 * kept * product(whole / length)
                + integrate(length - whole, |back| {
                    let s = 1.0 - back / length;
                    integrate(back / k, |u| conic_depth(1.0 + s, 0.5 + s, w, u))
                })
        };
        let exact = solid + added;
        let mut bits = Vec::new();
        for (face, a, b) in [(&top, &on_top, &on_wall), (&wall, &on_wall, &on_top)] {
            let built = ogeom_fillet::fillet_edge_distances_variable(
                &mut m, &s, &edge, face, a, b, profile, T,
            )
            .unwrap();
            let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{lean}: {:?}", diagnosis.problems);
            let v = volume(&m, &built.shape);
            assert!(
                (v - exact).abs() <= 1e-6,
                "{lean}: varying, volume {v} against {exact}"
            );
            let tangency = tangent_to_top_and_wall(&m, &built.shape);
            assert!(
                tangency <= 1e-7,
                "{lean}: the fill meets a face at {tangency}"
            );
            bits.push(v.to_bits());
        }
        assert_eq!(bits[0], bits[1], "{lean}: naming the other face first");
    }
}

/// Gauss-Legendre over `[lo, hi]` in 100 pieces, the variable laid as
/// a smoothstep so square-root ends at either side integrate smooth.
fn integrate_between(lo: f64, hi: f64, f: impl Fn(f64) -> f64) -> f64 {
    let nodes = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
    ];
    if hi <= lo {
        return 0.0;
    }
    let mut sum = 0.0;
    for k in 0..100 {
        let (a, b) = (f64::from(k) / 100.0, f64::from(k + 1) / 100.0);
        for (x, wt) in nodes {
            let s = (b - a).mul_add(x, a + b) * 0.5;
            let t = s * s * 2.0f64.mul_add(-s, 3.0);
            let jacobian = 6.0 * s * (1.0 - s);
            sum += wt * (b - a) * 0.5 * f((hi - lo).mul_add(t, lo)) * (hi - lo) * jacobian;
        }
    }
    sum
}

/// How far along x a block's top reaches at a given y.
type Reach = fn(f64) -> f64;

/// A 30 x 20 x 10 block whose end is a drum standing on the z axis, under
/// the 45 x 20 x 5 slab from (-5, 5, 10): the concave edge along
/// (x, 5, 10) from x = 0 to the drum. Leaning out, the block is cut to
/// within the drum of radius 20 about (5, 0); leaning in, the drum of
/// radius 12 about (35, 0) is cut out of it. With the edge, the top face,
/// the slab's front, how far along x the top reaches at `y` (the drum's
/// line there), and the solid's volume in closed form.
fn slab_on_drum_ended_block(m: &mut Model, out: bool) -> (Shape, Shape, Shape, Shape, Reach, f64) {
    let block = ogeom_algo::make_box(m, Frame::WORLD, (30.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let (centre, radius) = if out { (5.0, 20.0) } else { (35.0, 12.0) };
    let at = Frame::new(Point::new(centre, 0.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drum = ogeom_algo::make_cylinder(m, at, radius, 12.0, T)
        .unwrap()
        .shape;
    let base = if out {
        ogeom_bool::common(m, &block, &drum, T).unwrap().shape
    } else {
        ogeom_bool::cut(m, &block, &drum, T).unwrap().shape
    };
    let slab = ogeom_algo::make_box(
        m,
        Frame::new(Point::new(-5.0, 5.0, 10.0), Direction::Z, Direction::X, T).unwrap(),
        (45.0, 20.0, 5.0),
        T,
    )
    .unwrap()
    .shape;
    let s = ogeom_bool::fuse(m, &base, &slab, T).unwrap().shape;
    let edge = edge_where(m, &s, |p| {
        (p.y - 5.0).abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9
    });
    let top = face_at(m, &s, Point::new(10.0, 2.5, 10.0));
    let wall = face_at(m, &s, Point::new(10.0, 5.0, 12.5));
    // The block's plan: within the drum, the area under x = 5 + sqrt(400 -
    // y^2) over 0 <= y <= 20 less the sliver where its left side passes
    // x = 0; outside the drum, the rectangle less the segment the drum
    // cuts off past x = 23.
    let primitive =
        |r: f64, t: f64| 0.5 * t * r.mul_add(r, -t * t).sqrt() + 0.5 * r * r * (t / r).asin();
    let plan = if out {
        let start = 375.0f64.sqrt();
        100.0 + primitive(20.0, 20.0)
            - (5.0 * (20.0 - start) - (primitive(20.0, 20.0) - primitive(20.0, start)))
    } else {
        600.0 - (primitive(12.0, -5.0) - primitive(12.0, -12.0))
    };
    let reach: Reach = if out {
        |y| 5.0 + y.mul_add(-y, 400.0).sqrt()
    } else {
        |y| 35.0 - y.mul_add(-y, 144.0).sqrt()
    };
    (s, edge, top, wall, reach, 10.0f64.mul_add(plan, 4500.0))
}

/// A concave edge ending at a drum: the fill ends in the drum, as it ends
/// in a plane not square to it. On the slab on a block whose end is a
/// drum (leaning out or in), each point of the fill's section `u` into
/// the top from the crease runs along x from the block's square end to
/// the drum's line at y = 5 - u, run on past the edge's end where the
/// drum leans out, the laws held level there. With constant distances
/// that is the depth at `u` times that length, integrated over the top's
/// leg; with the laws 1 to 2 on the top and 0.5 to 1.5 on the slab (laid
/// from the edge's start), each point's depth under its station's laws
/// integrated over the plan. The same solid whichever face is named;
/// tangent to the top and the slab. An end face that is neither a plane
/// nor a cylinder is refused by name, and so are varying distances on a
/// D-shaped boss's straight foot, whose drum lies beyond the flat and so
/// does not span the fill's end.
#[test]
fn a_concave_edge_ending_at_a_drum_fills_up_to_the_drum() {
    use ogeom_fillet::RadiusLaw;
    use ogeom_geom::Curve3d as _;
    let rho = 0.4;
    let w = conic(rho);
    let profile = RoundProfile::Conic(rho);
    for out in [true, false] {
        let mut m = Model::new();
        let (s, edge, top, wall, reach, solid) = slab_on_drum_ended_block(&mut m, out);
        let base = volume(&m, &s);
        assert!(
            (base - solid).abs() <= 1e-8,
            "{out}: the block {base} against {solid}"
        );
        let length = reach(5.0);
        let (curve, (lo, _)) = ogeom_algo::placed_edge_curve(&m, &edge, T)
            .unwrap()
            .unwrap();
        let from_start = curve.point_at(lo, T).unwrap().x > 0.5 * length;
        let station = |x: f64| {
            let s = (x / length).clamp(0.0, 1.0);
            if from_start { 1.0 - s } else { s }
        };

        // Constant distances.
        let exact = solid + integrate(1.5, |u| conic_depth(1.5, 1.0, w, u) * reach(5.0 - u));
        let mut bits = Vec::new();
        for (face, a, b) in [(&top, 1.5, 1.0), (&wall, 1.0, 1.5)] {
            let built =
                ogeom_fillet::fillet_edge_distances(&mut m, &s, &edge, face, a, b, profile, T)
                    .unwrap();
            let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{out}: {:?}", diagnosis.problems);
            let v = volume(&m, &built.shape);
            assert!(
                (v - exact).abs() <= 1e-8,
                "{out}: constant, volume {v} against {exact}"
            );
            let tangency = tangent_to_top_and_wall(&m, &built.shape);
            assert!(
                tangency <= 1e-7,
                "{out}: the fill meets a face at {tangency}"
            );
            bits.push(v.to_bits());
        }
        assert_eq!(bits[0], bits[1], "{out}: naming the other face first");

        // Distances varying, each over its fraction of the edge from its
        // start.
        let on_top = RadiusLaw::Linear {
            start: 1.0,
            end: 2.0,
        };
        let on_wall = RadiusLaw::Linear {
            start: 0.5,
            end: 1.5,
        };
        let depth = |x: f64, u: f64| {
            let s = station(x);
            conic_depth(1.0 + s, 0.5 + s, w, u)
        };
        // Along x at depth u: the laws' ends, where the top's leg reaches
        // u, and the drum.
        let along = |u: f64| {
            let end = reach(5.0 - u);
            let mut cuts = vec![0.0, end];
            if end > length {
                cuts.push(length);
            }
            if (1.0..2.0).contains(&u) {
                let x = (u - 1.0) * length;
                cuts.push(if from_start { length - x } else { x });
            }
            cuts.retain(|x| (0.0..=end).contains(x));
            cuts.sort_by(f64::total_cmp);
            cuts.windows(2)
                .map(|p| integrate_between(p[0], p[1], |x| depth(x, u)))
                .sum::<f64>()
        };
        let exact = solid + integrate_between(0.0, 1.0, along) + integrate_between(1.0, 2.0, along);
        let mut bits = Vec::new();
        for (face, a, b) in [(&top, &on_top, &on_wall), (&wall, &on_wall, &on_top)] {
            let built = ogeom_fillet::fillet_edge_distances_variable(
                &mut m, &s, &edge, face, a, b, profile, T,
            )
            .unwrap();
            let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{out}: {:?}", diagnosis.problems);
            let v = volume(&m, &built.shape);
            assert!(
                (v - exact).abs() <= 1e-7,
                "{out}: varying, volume {v} against {exact}"
            );
            let tangency = tangent_to_top_and_wall(&m, &built.shape);
            assert!(
                tangency <= 1e-7,
                "{out}: the fill meets a face at {tangency}"
            );
            bits.push(v.to_bits());
        }
        assert_eq!(bits[0], bits[1], "{out}: naming the other face first");
    }

    // The block's end a ball instead: refused by name, constant or varying.
    let mut m = Model::new();
    let block = ogeom_algo::make_box(&mut m, Frame::WORLD, (30.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(35.0, 0.0, 5.0), Direction::Z, Direction::X, T).unwrap();
    let ball = ogeom_algo::make_sphere(&mut m, at, 12.0, T).unwrap().shape;
    let base = ogeom_bool::cut(&mut m, &block, &ball, T).unwrap().shape;
    let slab = ogeom_algo::make_box(
        &mut m,
        Frame::new(Point::new(-5.0, 5.0, 10.0), Direction::Z, Direction::X, T).unwrap(),
        (45.0, 20.0, 5.0),
        T,
    )
    .unwrap()
    .shape;
    let s = ogeom_bool::fuse(&mut m, &base, &slab, T).unwrap().shape;
    let edge = edge_where(&m, &s, |p| {
        (p.y - 5.0).abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9
    });
    let top = face_at(&m, &s, Point::new(10.0, 2.5, 10.0));
    let refused = |r: ogeom_core::OgeomResult<ogeom_algo::Built>| {
        let message = r.expect_err("refused").to_string();
        assert!(
            message.contains("neither a plane nor a cylinder"),
            "{message}"
        );
    };
    refused(ogeom_fillet::fillet_edge_distances(
        &mut m, &s, &edge, &top, 1.5, 1.0, profile, T,
    ));
    let law = |start, end| RadiusLaw::Linear { start, end };
    refused(ogeom_fillet::fillet_edge_distances_variable(
        &mut m,
        &s,
        &edge,
        &top,
        &law(1.0, 2.0),
        &law(0.5, 1.5),
        profile,
        T,
    ));

    let (mut m, boss, ..) = half_drum();
    let plate = ogeom_algo::make_box(
        &mut m,
        Frame::new(
            Point::new(-10.0, -10.0, -2.0),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap(),
        (20.0, 20.0, 2.0),
        T,
    )
    .unwrap()
    .shape;
    let s = ogeom_bool::fuse(&mut m, &plate, &boss, T).unwrap().shape;
    let foot = edge_where(&m, &s, |p| p.x.abs() < 1e-9 && p.z.abs() < 1e-9);
    let floor = face_at(&m, &s, Point::new(-5.0, 0.0, 0.0));
    let message = ogeom_fillet::fillet_edge_distances_variable(
        &mut m,
        &s,
        &foot,
        &floor,
        &law(1.0, 1.5),
        &law(0.5, 1.0),
        profile,
        T,
    )
    .expect_err("refused")
    .to_string();
    assert!(
        message.contains("does not span the fill's end"),
        "{message}"
    );
}

/// The largest angle a face of `shape` on a fitted surface meets a
/// neighbour at, sampled along each edge they share, but along the line of
/// the L-bracket's upright front edge (x = 8, y = 0): an edge left sharp
/// stays sharp up to the vertex.
fn tangent_but_at_the_sharp_edge(m: &Model, shape: &Shape) -> f64 {
    use ogeom_geom::Curve3d as _;
    let mut worst = 0.0_f64;
    for f in explore_unique(m, shape, ShapeType::Face).unwrap() {
        let data = m.node(&f).unwrap().data().as_face().unwrap();
        if !matches!(
            m.geometry().surface(data.surface),
            Some(ogeom_geom::SurfaceGeometry::BSpline(_))
        ) {
            continue;
        }
        for c in ogeom_fillet::analyse_blend(m, shape, &f, 17, T).unwrap() {
            let on_sharp = ogeom_algo::placed_edge_curve(m, &c.edge, T)
                .unwrap()
                .is_some_and(|(curve, (lo, hi))| {
                    [lo, f64::midpoint(lo, hi), hi].iter().all(|t| {
                        let p = curve.point_at(*t, T).unwrap();
                        (p.x - 8.0).abs() < 1e-6 && p.y.abs() < 1e-6
                    })
                });
            if !on_sharp {
                worst = worst.max(c.tangency_error);
            }
        }
    }
    worst
}

/// How far any point of `shape`'s mesh within 10 of the L-bracket's
/// front corner (8, 0, 8) stands outside the bracket and the fill of the
/// conic corner of legs `fill` (along the floor, up the wall) and weight
/// `w` on its re-entrant edge. A convex band only takes material off.
fn outside_bracket_and_fill(m: &Model, shape: &Shape, fill: [f64; 2], w: f64) -> f64 {
    let deflection = ogeom_mesh::Deflection {
        chord: 1e-3,
        ..ogeom_mesh::Deflection::default()
    };
    let out = |p: Point| -> f64 {
        let q = [p.x, p.y, p.z];
        let in_box = |lo: [f64; 3], hi: [f64; 3]| {
            (0..3)
                .map(|i| (lo[i] - q[i]).max(q[i] - hi[i]))
                .fold(f64::NEG_INFINITY, f64::max)
        };
        let floor = in_box([0.0, 0.0, 0.0], [20.0, 20.0, 8.0]);
        let upright = in_box([0.0, 0.0, 0.0], [8.0, 20.0, 20.0]);
        let (u, v) = (p.x - 8.0, p.z - 8.0);
        let under = if u < fill[0] && p.y >= -1e-9 {
            v - conic_depth(fill[0], fill[1], w, u.max(0.0))
        } else {
            f64::INFINITY
        };
        floor.min(upright).min(under)
    };
    let mut worst = f64::NEG_INFINITY;
    for f in explore_unique(m, shape, ShapeType::Face).unwrap() {
        for p in ogeom_mesh::triangulate_face(m, &f, deflection, T)
            .unwrap()
            .positions
        {
            if p.distance(Point::new(8.0, 0.0, 8.0)) <= 10.0 {
                worst = worst.max(out(p));
            }
        }
    }
    worst
}

/// The L-bracket's front corner (8, 0, 8): its floor's front edge convex,
/// its re-entrant edge concave and the upright's front edge left sharp.
/// With an asymmetric round on either or both of the first two, the
/// convex band cuts and the concave one fills, and one patch closes the
/// corner between them, tangent to both bands but along the sharp edge,
/// which stays sharp up to the vertex. The same solid in either order and
/// whichever face is named. Away from the corner the bracket loses the
/// convex corner along the floor's 12 and gains the fill along the 20 of
/// the re-entrant edge; what the corner itself adds or takes stays within
/// the cube of the widest band's widest leg. Every point near the corner
/// lies in the bracket or the fill.
#[test]
#[ignore = "heavy: nine corner patches and a fine mesh of each"]
fn a_convex_and_a_concave_band_beside_an_edge_left_sharp_close_on_one_patch() {
    let w = conic(0.4);
    let asym = |face: &Shape, a: f64, b: f64| Blend::RoundDistances {
        face: face.clone(),
        on_face: a,
        on_other: b,
        profile: RoundProfile::Conic(0.4),
    };
    let mut m = Model::new();
    let (s, crease, floor, wall) = l_bracket(&mut m);
    let front_edge = edge_where(&m, &s, |p| {
        p.y.abs() < 1e-9 && (p.z - 8.0).abs() < 1e-9 && p.x >= 8.0 - 1e-9
    });
    let front = face_at(&m, &s, Point::new(14.0, 0.0, 4.0));
    // A circular round of radius r takes r^2 (1 - pi / 4) off a square
    // corner, or adds it, its legs r each way.
    let round = |r: f64| r * r * (1.0 - core::f64::consts::FRAC_PI_4);
    /// One request at the corner: each band named one way and the other,
    /// what the convex band takes and the concave one adds per unit length,
    /// the fill's legs (along the floor, up the wall) and its weight, and
    /// the widest leg.
    struct Case {
        convex: [Blend; 2],
        concave: [Blend; 2],
        taken: f64,
        filled: f64,
        fill: ([f64; 2], f64),
        widest: f64,
    }
    let cases = [
        Case {
            convex: [asym(&floor, 2.0, 1.0), asym(&front, 1.0, 2.0)],
            concave: [Blend::Round(1.5), Blend::Round(1.5)],
            taken: corner_area(2.0, 1.0, FRAC_PI_2, w),
            filled: round(1.5),
            fill: ([1.5, 1.5], elliptical(FRAC_PI_2)),
            widest: 2.0,
        },
        Case {
            convex: [Blend::Round(1.5), Blend::Round(1.5)],
            concave: [asym(&floor, 2.0, 1.0), asym(&wall, 1.0, 2.0)],
            taken: round(1.5),
            filled: corner_area(2.0, 1.0, FRAC_PI_2, w),
            fill: ([2.0, 1.0], w),
            widest: 2.0,
        },
        Case {
            convex: [asym(&floor, 1.2, 1.8), asym(&front, 1.8, 1.2)],
            concave: [asym(&floor, 1.6, 1.0), asym(&wall, 1.0, 1.6)],
            taken: corner_area(1.2, 1.8, FRAC_PI_2, w),
            filled: corner_area(1.6, 1.0, FRAC_PI_2, w),
            fill: ([1.6, 1.0], w),
            widest: 1.8,
        },
    ];
    for (k, case) in cases.into_iter().enumerate() {
        let Case {
            convex: [convex, convex_other],
            concave: [concave, concave_other],
            taken,
            filled,
            fill: (fill, weight),
            widest,
        } = case;
        let request = [(front_edge.clone(), convex), (crease.clone(), concave)];
        let v = either_order(&mut m, &s, &request);
        let other = [
            (front_edge.clone(), convex_other),
            (crease.clone(), concave_other),
        ];
        let again = ogeom_fillet::blend_edges(&mut m, &s, &other, ChamferCorner::Point, T).unwrap();
        assert_eq!(
            volume(&m, &again.shape).to_bits(),
            v.to_bits(),
            "case {k}: naming the other faces"
        );
        let built =
            ogeom_fillet::blend_edges(&mut m, &s, &request, ChamferCorner::Point, T).unwrap();
        let tangency = tangent_but_at_the_sharp_edge(&m, &built.shape);
        assert!(
            tangency <= 2e-4,
            "case {k}: the patch meets a band at {tangency} rad"
        );
        let past = outside_bracket_and_fill(&m, &built.shape, fill, weight);
        assert!(
            past <= 2e-4,
            "case {k}: the corner stands {past} outside the bracket and its fill"
        );
        let near = v - (5120.0 - 12.0 * taken + 20.0 * filled);
        let cube = widest.powi(3);
        assert!(
            near.abs() <= cube,
            "case {k}: {near} at the corner, past {cube}"
        );
    }
}

/// The blend a bevel names on `face`: `a` along it, `b` along the other.
fn bevel(face: &Shape, a: f64, b: f64) -> Blend {
    Blend::Chamfer(ogeom_fillet::Chamfer::Distances {
        face: face.clone(),
        on_face: a,
        on_other: b,
    })
}

/// A bevel's depth `u` from its edge along the face of its leg `leg`, its
/// other leg `height`: falling straight to nothing at the contact.
fn falling(leg: f64, height: f64, u: f64) -> f64 {
    (height * (1.0 - u / leg)).max(0.0)
}

/// A bevel and an asymmetric round of one sense at a vertex whose third
/// edge is left sharp meet by intersection, as two rounds do: on the box's
/// top corner the round of the top's front edge and the bevel of its left
/// edge each cut along their whole edge, and each point of the top under
/// both takes the shallower depth. The same solid in either order.
#[test]
fn a_bevel_and_an_asymmetric_round_beside_an_edge_left_sharp_meet_by_intersection() {
    let w = conic(0.4);
    let mut m = Model::new();
    let (s, [x, y, _], [top, ..]) = box_corner(&mut m);
    let (a_top, b_front) = (3.0, 1.5);
    let (b_top, a_left) = (1.0, 2.0);
    let both = integrate(b_top, |u| {
        let k = falling(b_top, a_left, u);
        integrate(a_top, |t| conic_depth(a_top, b_front, w, t).min(k))
    });
    let exact =
        4000.0 - 20.0 * corner_area(a_top, b_front, FRAC_PI_2, w) - 20.0 * 0.5 * a_left * b_top
            + both;
    let round = Blend::RoundDistances {
        face: top.clone(),
        on_face: a_top,
        on_other: b_front,
        profile: RoundProfile::Conic(0.4),
    };
    let v = either_order(&mut m, &s, &[(x, round), (y, bevel(&top, b_top, a_left))]);
    assert!((v - exact).abs() <= 1e-5, "volume {v} against {exact}");
}

/// At the pocket's inside corner, the upright left sharp, an asymmetric
/// round's fill and a bevel's each add along their edge, less once the
/// shallower over the floor under both. The same solid in either order.
#[test]
fn a_bevel_and_an_asymmetric_fill_beside_an_edge_left_sharp_meet_by_intersection() {
    let w = conic(0.4);
    let mut m = Model::new();
    let (s, [ex, ey, _], [floor, ..]) = pocket(&mut m);
    let (a_floor, b_wall) = (1.5, 1.0);
    let (c_floor, d_wall) = (0.8, 1.2);
    let both = integrate(c_floor, |u| {
        let k = falling(c_floor, d_wall, u);
        integrate(a_floor, |t| conic_depth(a_floor, b_wall, w, t).min(k))
    });
    let exact =
        3500.0 + 10.0 * corner_area(a_floor, b_wall, FRAC_PI_2, w) + 10.0 * 0.5 * c_floor * d_wall
            - both;
    let fill = Blend::RoundDistances {
        face: floor.clone(),
        on_face: a_floor,
        on_other: b_wall,
        profile: RoundProfile::Conic(0.4),
    };
    let v = either_order(
        &mut m,
        &s,
        &[(ex, fill), (ey, bevel(&floor, c_floor, d_wall))],
    );
    assert!((v - exact).abs() <= 1e-5, "volume {v} against {exact}");
}

/// A bevel and an asymmetric round of the other sense at a vertex whose
/// third edge is left sharp: bevels first (ADR 0001, section 13,
/// question 3), the round then on what the bevel left, where its edge ends
/// in the bevel's plane. On the L-bracket, a bevel of the floor's front
/// edge `p` along the floor and `q` down the front, then a fill of the
/// crease: each point of the fill's section `h` above the floor stands
/// behind the bevel's plane only from `p + p h / q` along the crease, so
/// the fill adds `(20 - p) A - (p / q) M`, `M` the section's first moment
/// about the floor. A bevel filling the crease `p` along the floor and `q`
/// up the wall, then a round of the floor's front edge: the band, cut by
/// the bevel's plane carried down into the floor, takes each point `t`
/// below the floor from `p + p t / q` past the wall, `(12 - p) A - (p / q)
/// M`. The same solid in either order and naming the other faces.
#[test]
fn an_asymmetric_round_and_a_bevel_of_the_other_sense_beside_an_edge_left_sharp_take_the_bevel_first()
 {
    let w = conic(0.4);
    let asym = |face: &Shape, a: f64, b: f64| Blend::RoundDistances {
        face: face.clone(),
        on_face: a,
        on_other: b,
        profile: RoundProfile::Conic(0.4),
    };
    let (a, b) = (1.5, 1.0);
    let (p, q) = (1.2, 0.8);
    let area = corner_area(a, b, FRAC_PI_2, w);
    let moment = integrate(a, |u| 0.5 * conic_depth(a, b, w, u).powi(2));
    let mut m = Model::new();
    let (s, crease, floor, wall) = l_bracket(&mut m);
    let front_edge = edge_where(&m, &s, |p| {
        p.y.abs() < 1e-9 && (p.z - 8.0).abs() < 1e-9 && p.x >= 8.0 - 1e-9
    });
    let front = face_at(&m, &s, Point::new(14.0, 0.0, 4.0));
    let cases = [
        (
            [
                (front_edge.clone(), bevel(&floor, p, q)),
                (crease.clone(), asym(&floor, a, b)),
            ],
            [
                (front_edge.clone(), bevel(&front, q, p)),
                (crease.clone(), asym(&wall, b, a)),
            ],
            5120.0 - 12.0 * 0.5 * p * q + (20.0 - p) * area - p / q * moment,
        ),
        (
            [
                (crease.clone(), bevel(&floor, p, q)),
                (front_edge.clone(), asym(&floor, a, b)),
            ],
            [
                (crease.clone(), bevel(&wall, q, p)),
                (front_edge.clone(), asym(&front, b, a)),
            ],
            5120.0 + 20.0 * 0.5 * p * q - ((12.0 - p) * area - p / q * moment),
        ),
    ];
    for (k, (request, other, exact)) in cases.into_iter().enumerate() {
        let v = either_order(&mut m, &s, &request);
        assert!(
            (v - exact).abs() <= 1e-7,
            "case {k}: volume {v} against {exact}"
        );
        let again = ogeom_fillet::blend_edges(&mut m, &s, &other, ChamferCorner::Point, T).unwrap();
        assert_eq!(
            volume(&m, &again.shape).to_bits(),
            v.to_bits(),
            "case {k}: naming the other faces"
        );
    }
}

/// A 20 x 12 x 4 block, its upright at (20, 0) rounded at 3: the top's
/// front edge to x = 17, the quarter arc round (17, 3) and the top's right
/// edge from y = 3, tangent end to end, with the top face. The chain in
/// that order and its length.
fn rounded_corner_chain(m: &mut Model) -> (Shape, [Shape; 3], Shape, f64) {
    let block = ogeom_algo::make_box(m, Frame::WORLD, (20.0, 12.0, 4.0), T)
        .unwrap()
        .shape;
    let upright = edge_where(m, &block, |p| (p.x - 20.0).abs() < 1e-9 && p.y.abs() < 1e-9);
    let s = ogeom_fillet::fillet_edge(m, &block, &upright, 3.0, T)
        .unwrap()
        .shape;
    let top = face_at(m, &s, Point::new(10.0, 6.0, 4.0));
    let front = edge_where(m, &s, |p| p.y.abs() < 1e-9 && (p.z - 4.0).abs() < 1e-9);
    let right = edge_where(m, &s, |p| {
        (p.x - 20.0).abs() < 1e-9 && (p.z - 4.0).abs() < 1e-9
    });
    let arc = edge_where(m, &s, |p| {
        (p.z - 4.0).abs() < 1e-9
            && ((p.x - 17.0).abs() < 1e-9 && p.y.abs() < 1e-9
                || (p.x - 20.0).abs() < 1e-9 && (p.y - 3.0).abs() < 1e-9)
    });
    let length = 17.0 + 1.5 * core::f64::consts::PI + 9.0;
    (s, [front, arc, right], top, length)
}

/// The tangent planes' gap at a junction of the rounded corner's chain
/// `sigma` along it, of the round whose sections stand square to the
/// crease: legs `d(sigma)` (along the top, down the wall), weight `w`.
/// Where the laws vary, a section's points move along the crease at a pace
/// that jumps at the junction by how far each point stands in from the
/// crease over the arc's radius, while the laws move them across it alike
/// on both sides; so the two bands' tangent planes part there by an angle
/// of the laws' slope times that jump. The largest over the section, each
/// side's along-derivative taken one-sidedly.
fn junction_gap(sigma: f64, d: impl Fn(f64) -> (f64, f64), w: f64) -> f64 {
    use core::f64::consts::PI;
    let arc = 1.5 * PI;
    // The crease point and the top's leg inward, by arc length.
    let crease = |s: f64| -> (Point, Vector) {
        if s <= 17.0 {
            (Point::new(s, 0.0, 4.0), Y)
        } else if s <= 17.0 + arc {
            let th = (s - 17.0) / 3.0;
            (
                Point::new(17.0 + 3.0 * th.sin(), 3.0 - 3.0 * th.cos(), 4.0),
                Vector::new(-th.sin(), th.cos(), 0.0),
            )
        } else {
            (Point::new(20.0, 3.0 + s - 17.0 - arc, 4.0), -X)
        }
    };
    let section = |s: f64, t: f64| -> Vector {
        let (c, u) = crease(s);
        let (a, b) = d(s);
        let c = c - Point::ORIGIN;
        let (p0, p2) = (c + u * a, c - Vector::new(0.0, 0.0, b));
        let (b0, b1, b2) = ((1.0 - t) * (1.0 - t), 2.0 * t * (1.0 - t) * w, t * t);
        (p0 * b0 + c * b1 + p2 * b2) * (1.0 / (b0 + b1 + b2))
    };
    let h = 1e-6;
    let mut worst = 0.0_f64;
    for k in 1..64 {
        let t = f64::from(k) / 64.0;
        let p = section(sigma, t);
        let across = section(sigma, t + h) - section(sigma, t - h);
        let n = |v: Vector| v * (1.0 / v.magnitude());
        let before = n(across.cross(p - section(sigma - h, t)));
        let after = n(across.cross(section(sigma + h, t) - p));
        worst = worst.max(before.dot(after).clamp(-1.0, 1.0).acos());
    }
    worst
}

/// Each place two bands of `s` meet (a band being a face on neither a
/// plane nor a cylinder), as a point of the edge between them, and the
/// angle between their normals along it.
fn band_joins(m: &Model, s: &Shape) -> Vec<(Point, f64)> {
    use ogeom_geom::Curve3d as _;
    let band = |f: &Shape| {
        let data = m.node(f).unwrap().data().as_face().unwrap();
        !matches!(
            m.geometry().surface(data.surface),
            Some(ogeom_geom::SurfaceGeometry::Plane(_) | ogeom_geom::SurfaceGeometry::Cylinder(_))
        )
    };
    let mut joins = Vec::new();
    for f in explore_unique(m, s, ShapeType::Face).unwrap() {
        if !band(&f) {
            continue;
        }
        for c in ogeom_fillet::analyse_blend(m, s, &f, 65, T).unwrap() {
            if band(&c.neighbour) {
                let (curve, (lo, hi)) = ogeom_algo::placed_edge_curve(m, &c.edge, T)
                    .unwrap()
                    .unwrap();
                let p = curve.point_at(f64::midpoint(lo, hi), T).unwrap();
                joins.push((p, c.tangency_error));
            }
        }
    }
    joins
}

/// Each band-to-band join of a round along the rounded corner's chain held
/// to the model's own gap at its junction, within `within`; both junctions
/// met, from either side.
fn junctions_hold(
    m: &Model,
    s: &Shape,
    d: impl Fn(f64) -> (f64, f64) + Copy,
    w: f64,
    within: f64,
    label: &str,
) {
    use core::f64::consts::PI;
    let at = [17.0, 17.0 + 1.5 * PI];
    let joins = band_joins(m, s);
    let mut met = [0; 2];
    for (p, measured) in joins {
        // The first junction's section stands in x = 17, the second's in
        // y = 3.
        let k = usize::from(p.x > 17.5 && (p.y - 3.0).abs() < 1e-6);
        assert!(
            k == 1 || (p.x - 17.0).abs() < 1e-6,
            "{label}: bands meet at {p:?}, at no junction"
        );
        met[k] += 1;
        let model = junction_gap(at[k], d, w);
        assert!(
            (measured - model).abs() <= within,
            "{label}: the bands at {p:?} meet at {measured:.3e} rad, the model's sections at \
             {model:.3e}"
        );
    }
    assert!(met == [2, 2], "{label}: joins {met:?} at the two junctions");
}

/// The chain round the rounded corner asked with one law on the top, 0.6
/// to 1.2, and one on the walls, 0.4 to 1.0, laid along the whole chain by
/// arc length. Each station's section is the square conic corner of its
/// distances: along the straight edges it takes that corner's area,
/// round the arc its area times `(R - c) / R`, `c` its centroid's depth
/// into the top, by Pappus station by station; the laws run from the
/// front edge's start the way it runs. Valid, the edges consumed, each
/// band tangent to the faces it meets within the march's fit, and within
/// that fit of what the band takes. At each junction the bands meet as
/// the sections square to the crease do ([`junction_gap`]), within the
/// march's fit. The walls are three faces, so naming the front alone is
/// refused by name.
#[test]
fn distances_varying_along_a_tangent_chain_take_the_integral_of_their_corners() {
    use core::f64::consts::PI;
    use ogeom_fillet::RadiusLaw;
    let on_top = RadiusLaw::Linear {
        start: 0.6,
        end: 1.2,
    };
    let on_wall = RadiusLaw::Linear {
        start: 0.4,
        end: 1.0,
    };
    for profile in PROFILES {
        let w = weight_of(profile);
        let mut m = Model::new();
        let (s, chain, top, length) = rounded_corner_chain(&mut m);
        let before = volume(&m, &s);
        // Where the laws start: the front edge's first vertex as it runs.
        let (first, _) = ogeom_algo::edge_vertices(&m, &chain[0]).unwrap().unwrap();
        let from_left = point_of(&m, &first).x < 1.0;
        let at = |sigma: f64| {
            let t = sigma / length;
            let s = if from_left { t } else { 1.0 - t };
            (0.6 + 0.6 * s, 0.4 + 0.6 * s)
        };
        let arc = 1.5 * PI;
        let straight = |from: f64, run: f64| {
            integrate(run, |u| {
                let (a, b) = at(from + u);
                corner_area(a, b, FRAC_PI_2, w)
            })
        };
        let turned = integrate(arc, |u| {
            let (a, b) = at(17.0 + u);
            let (area, c) = square_corner(a, b, w);
            area * (3.0 - c[0]) / 3.0
        });
        let taken = straight(0.0, 17.0) + turned + straight(17.0 + arc, 9.0);
        let exact = before - taken;
        let built = ogeom_fillet::fillet_edges_distances_variable(
            &mut m, &s, &chain, &top, &on_top, &on_wall, profile, T,
        )
        .unwrap();
        let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "{profile:?}: {:?}",
            diagnosis.problems
        );
        for edge in &chain {
            assert!(built.history.is_deleted(edge));
            assert!(!built.history.generated(edge).is_empty(), "each band");
        }
        let v = volume(&m, &built.shape);
        assert!(
            (v - exact).abs() <= 1e-4 * taken,
            "{profile:?}: volume {v} against {exact}"
        );
        let mut worst = 0.0_f64;
        for f in explore_unique(&m, &built.shape, ShapeType::Face).unwrap() {
            let data = m.node(&f).unwrap().data().as_face().unwrap();
            if matches!(
                m.geometry().surface(data.surface),
                Some(
                    ogeom_geom::SurfaceGeometry::Plane(_)
                        | ogeom_geom::SurfaceGeometry::Cylinder(_)
                )
            ) {
                continue;
            }
            for c in ogeom_fillet::analyse_blend(&m, &built.shape, &f, 17, T).unwrap() {
                use ogeom_geom::Curve3d as _;
                let (curve, (lo, hi)) = ogeom_algo::placed_edge_curve(&m, &c.edge, T)
                    .unwrap()
                    .unwrap();
                let p = curve.point_at(f64::midpoint(lo, hi), T).unwrap();
                // Along the faces it meets, not the block's end faces its
                // ends stand in, square to it.
                let nd = m.node(&c.neighbour).unwrap().data().as_face().unwrap();
                let on_a_face = matches!(
                    m.geometry().surface(nd.surface),
                    Some(
                        ogeom_geom::SurfaceGeometry::Plane(_)
                            | ogeom_geom::SurfaceGeometry::Cylinder(_)
                    )
                );
                if on_a_face && p.x.abs() > 1e-7 && (p.y - 12.0).abs() > 1e-7 {
                    worst = worst.max(c.tangency_error);
                }
            }
        }
        assert!(
            worst <= 1e-3,
            "{profile:?}: a band meets a face at {worst} rad"
        );
        junctions_hold(&m, &built.shape, at, w, 1e-3, &format!("{profile:?}"));
        // Naming the walls: each straight edge's wall and the drum are
        // different faces, so the chain names the top only.
        let front = face_at(&m, &s, Point::new(10.0, 0.0, 2.0));
        let refused = ogeom_fillet::fillet_edges_distances_variable(
            &mut m, &s, &chain, &front, &on_wall, &on_top, profile, T,
        )
        .expect_err("refused")
        .to_string();
        assert!(refused.contains("beside every edge"), "{refused}");
    }
}

/// The same chain round the rounded corner under one radius law, 0.6 to
/// 1.2, laid along it by arc length: the circular round, each section the
/// square corner's quarter circle, its volume taken by Pappus as the
/// asymmetric chain's is. At each junction the bands meet as the sections
/// square to the crease do ([`junction_gap`]): the straight band exact
/// and the marched one within its fit of the model.
#[test]
fn a_radius_law_along_a_tangent_chain_meets_at_its_junctions_as_its_sections_do() {
    use core::f64::consts::PI;
    use ogeom_fillet::RadiusLaw;
    let mut m = Model::new();
    let (s, chain, _, length) = rounded_corner_chain(&mut m);
    let before = volume(&m, &s);
    let (first, _) = ogeom_algo::edge_vertices(&m, &chain[0]).unwrap().unwrap();
    let from_left = point_of(&m, &first).x < 1.0;
    let at = |sigma: f64| {
        let t = sigma / length;
        let r = 0.6 + 0.6 * if from_left { t } else { 1.0 - t };
        (r, r)
    };
    let w = FRAC_1_SQRT_2;
    let straight = |from: f64, run: f64| {
        integrate(run, |u| {
            let (a, b) = at(from + u);
            corner_area(a, b, FRAC_PI_2, w)
        })
    };
    let arc = 1.5 * PI;
    let turned = integrate(arc, |u| {
        let (a, b) = at(17.0 + u);
        let (area, c) = square_corner(a, b, w);
        area * (3.0 - c[0]) / 3.0
    });
    let taken = straight(0.0, 17.0) + turned + straight(17.0 + arc, 9.0);
    let law = RadiusLaw::Linear {
        start: 0.6,
        end: 1.2,
    };
    let built = ogeom_fillet::fillet_edges_variable(&mut m, &s, &chain, &law, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let v = volume(&m, &built.shape);
    assert!(
        (v - (before - taken)).abs() <= 1e-4 * taken,
        "volume {v} against {}",
        before - taken
    );
    junctions_hold(&m, &built.shape, at, w, 1e-5, "circular");
}

/// The drum standing in the slab's corner, its concave ruling at (5, 0, z)
/// asked with 0.2 along the slab's side per unit of scale and 1 round the
/// drum. The drum's leg turns `beta = s / 5`, and the drum's tangent there
/// meets the side's line `5 (1 / cos(beta) - 1)` short of the crease,
/// between the contacts: the shoulder the conic stands on reaches the
/// side's contact, `0.2 s` from the crease, long before the side's leg
/// leaves its face (at 25). So the largest scale is the root of
/// `5 (1 / cos(s / 5) - 1) = 0.2 s`, bisected here. The answer builds
/// valid, and a scale a hundredth past it is refused.
#[test]
fn the_largest_pair_on_a_concave_ruling_stops_where_the_shoulder_reaches_a_contact() {
    let mut m = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, 5.0, 4.0, T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(0.0, -10.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let slab = ogeom_algo::make_box(&mut m, at, (10.0, 10.0, 4.0), T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &drum, &slab, T).unwrap().shape;
    let crease = edge_where(&m, &s, |p| (p.x - 5.0).abs() < 1e-9 && p.y.abs() < 1e-9);
    let side = face_at(&m, &s, Point::new(8.0, 0.0, 2.0));
    let (on_side, on_drum) = (0.2, 1.0);
    let short = |s: f64| 5.0 * (1.0 / (on_drum * s / 5.0).cos() - 1.0) - on_side * s;
    let (mut lo, mut hi) = (0.5, 7.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if short(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let want = 0.5 * (lo + hi);
    let largest =
        ogeom_fillet::largest_fillet_distances(&m, &s, &crease, &side, on_side, on_drum, 100.0, T)
            .unwrap()
            .expect("a pair fits");
    assert!(
        (largest - want).abs() <= 1e-4 * want,
        "the largest scale {largest}, not {want}"
    );
    let profile = RoundProfile::Elliptical;
    let built = ogeom_fillet::fillet_edge_distances(
        &mut m,
        &s,
        &crease,
        &side,
        largest * on_side,
        largest * on_drum,
        profile,
        T,
    )
    .unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let past = ogeom_fillet::fillet_edge_distances(
        &mut m,
        &s,
        &crease,
        &side,
        largest * on_side * 1.01,
        largest * on_drum * 1.01,
        profile,
        T,
    );
    assert!(past.is_err(), "a scale past the largest is not refused");
}
