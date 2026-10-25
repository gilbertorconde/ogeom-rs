//! Pipe shells round sharp corners beyond straight legs under the
//! rotation-minimizing frame: corners against curved legs, closed spines,
//! and the binormal and auxiliary laws, measured against closed forms.
//!
//! A square section of side `2 a` with two sides along `z` swept along a
//! spine in the `xy` plane is the extrusion along `z`, through `2 a`, of
//! the plane region the spine's band of half width `a` covers; so each
//! volume is `2 a` times a plane area.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, edge_vertices, make_edge, make_edge_between, make_face, make_polygon, make_vertex,
    make_wire, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve3d as _, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::offset::{PipeCorners, PipeLaw, make_pipe_shell_with};
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// A square of side `2 a` about `at`, square to `heading` (in the `xy`
/// plane), two of its sides along `z`.
fn square_section(model: &mut Model, at: Point, heading: Vector, a: f64) -> Shape {
    rectangle_section(model, at, heading, a, a)
}

/// A rectangle `2 a` across and `2 b` along `z` about `at`, square to
/// `heading` (in the `xy` plane).
fn rectangle_section(model: &mut Model, at: Point, heading: Vector, a: f64, b: f64) -> Shape {
    let across = Vector::new(0.0, 0.0, 1.0).cross(heading);
    let across = across / across.magnitude();
    let up = Vector::new(0.0, 0.0, 1.0);
    let pts = [(-a, -b), (a, -b), (a, b), (-a, b)].map(|(s, t)| at + across * s + up * t);
    let wire = make_polygon(model, &pts, true, T).unwrap().shape;
    let frame = Frame::new(
        at,
        Direction::new(heading, T).unwrap(),
        Direction::new(across, T).unwrap(),
        T,
    )
    .unwrap();
    make_face(
        model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape
}

/// The volume, and whether it was integrated exactly.
fn measure(model: &Model, shape: &Shape) -> (f64, bool) {
    let p = volume_properties(model, shape, Deflection::with_chord(1e-4).unwrap(), T).unwrap();
    (p.mass, p.deflection == 0.0)
}

fn planes_only(model: &Model, shape: &Shape) -> bool {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .all(|f| {
            let Some(NodeData::Face(d)) = model.node(f).map(|n| n.data()) else {
                return false;
            };
            matches!(
                model.geometry().surface(d.surface),
                Some(SurfaceGeometry::Plane(_))
            )
        })
}

fn tilted() -> PipeLaw<'static> {
    PipeLaw::Binormal(Direction::new(Vector::new(1.0, 1.0, 1.0), T).unwrap())
}

/// A rectangle down the L (0,0,0) -> (0,L,0) -> (L,L,0) under a binormal
/// tilted off the L's plane: each leg's section is set by the law, not
/// turned across the corner, and the two legs meet on the mitre plane.
/// Each leg is a prism of its section cut square at its start and by the
/// mitre plane through the spine at the corner, which holds the section's
/// area times the leg's length along the spine; the legs lie either side
/// of the mitre plane.
#[test]
fn a_binormal_rectangle_down_an_l_mitres_to_its_legs_volume() {
    let (a, l) = (1.0, 10.0);
    let mut model = Model::new();
    let profile = rectangle_section(
        &mut model,
        Point::ORIGIN,
        Vector::new(0.0, 1.0, 0.0),
        a,
        2.0 * a,
    );
    let spine = make_polygon(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(0.0, l, 0.0),
            Point::new(l, l, 0.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let pipe = make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &tilted(),
        PipeCorners::Mitre,
        1e-3,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &pipe, T).unwrap();
    assert!(found.is_valid(), "{found}");
    assert!(planes_only(&model, &pipe));
    let (v, exact) = measure(&model, &pipe);
    assert!(exact);
    let want = 8.0 * a * a * 2.0 * l;
    assert!((v - want).abs() < want * 1e-12, "{v} against {want}");
}

/// The square ring (0,0) -> (L,0) -> (L,L) -> (0,L) -> (0,0), under the
/// tilted binormal and mitred at all four corners: the section's area
/// times the ring's length.
#[test]
fn a_binormal_rectangle_round_a_square_ring_mitres_to_its_length() {
    let (a, l) = (1.0, 10.0);
    let mut model = Model::new();
    let profile = rectangle_section(
        &mut model,
        Point::ORIGIN,
        Vector::new(1.0, 0.0, 0.0),
        a,
        2.0 * a,
    );
    let spine = make_polygon(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(l, 0.0, 0.0),
            Point::new(l, l, 0.0),
            Point::new(0.0, l, 0.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    let pipe = make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &tilted(),
        PipeCorners::Mitre,
        1e-3,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &pipe, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let (v, exact) = measure(&model, &pipe);
    assert!(exact);
    let want = 8.0 * a * a * 4.0 * l;
    assert!((v - want).abs() < want * 1e-12, "{v} against {want}");
}

/// An auxiliary spine running above the L keeps the section's normal
/// pointing up the whole way: the mitred L, its legs' volume.
#[test]
fn an_auxiliary_spine_above_an_l_mitres_to_its_legs_volume() {
    let (a, l, h) = (1.0, 10.0, 5.0);
    let mut model = Model::new();
    let profile = square_section(&mut model, Point::ORIGIN, Vector::new(0.0, 1.0, 0.0), a);
    let spine = make_polygon(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(0.0, l, 0.0),
            Point::new(l, l, 0.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let guide = make_polygon(
        &mut model,
        &[
            Point::new(0.0, -1.0, h),
            Point::new(0.0, l, h),
            Point::new(l + 1.0, l, h),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let pipe = make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::Auxiliary { guide: &guide },
        PipeCorners::Mitre,
        1e-3,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &pipe, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let (v, _) = measure(&model, &pipe);
    let want = 4.0 * a * a * 2.0 * l;
    assert!((v - want).abs() < want * 1e-9, "{v} against {want}");
}

/// The spine: a quarter arc of radius `r` about (0, r) from (-r, r),
/// heading -y, to the origin, heading +x; then a sharp quarter turn left
/// up the line to (0, l2).
fn arc_then_line(model: &mut Model, r: f64, l2: f64) -> Shape {
    let frame = Frame::new(
        Point::new(0.0, r, 0.0),
        Direction::Z,
        Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
        T,
    )
    .unwrap();
    let circle = Circle::new(frame, r, T).unwrap();
    let arc = make_edge(
        model,
        CircleCurve::new(circle).into(),
        (1.5 * PI, 2.0 * PI),
        T,
    )
    .unwrap()
    .shape;
    let (_, joint) = edge_vertices(model, &arc).unwrap().unwrap();
    let top = make_vertex(model, Point::new(0.0, l2, 0.0)).shape;
    let segment = LineCurve::segment(Point::ORIGIN, Point::new(0.0, l2, 0.0), T).unwrap();
    let range = segment.domain();
    let line_edge = make_edge_between(model, segment.into(), range, &joint, &top, T)
        .unwrap()
        .shape;
    make_wire(model, &[arc, line_edge], T).unwrap().shape
}

/// `∫₀ᵃ √(ρ² - x²) dx`.
fn under_circle(rho: f64, a: f64) -> f64 {
    a / 2.0 * (rho * rho - a * a).sqrt() + rho * rho / 2.0 * (a / rho).asin()
}

/// The arc's band is the annular sector between radii `r - a` and `r + a`,
/// of area `π r a`; the line's band from the corner is `2 a l2`.
///
/// Round: the line's band overlaps the arc's in `-a < x < 0` from `y = 0`
/// up to the inner circle, `a r - ∫₀ᵃ √((r - a)² - x²) dx`; the outer half
/// of the section turned about the corner adds a quarter disc of radius
/// `a` in `x > 0, y < 0`, which neither band reaches.
///
/// Extended: the arc's run on past the corner covers `0 < x < a` across
/// the band and the line's run back covers `-a < y < 0`, so the square
/// `[-a, a]²` is covered whole; past it the line's band adds
/// `2 a (l2 - a)`, and the arc's band adds all but its part over
/// `-a < x < 0`, `∫₀ᵃ (√((r + a)² - x²) - √((r - a)² - x²)) dx`.
#[test]
fn a_square_round_a_corner_after_an_arc_turns_each_way_to_its_closed_form() {
    let (a, r, l2) = (1.0, 8.0, 10.0);
    let round =
        2.0 * a * (PI * r * a + 2.0 * a * l2 + PI * a * a / 4.0 - (a * r - under_circle(r - a, a)));
    let extended = 2.0
        * a
        * (4.0 * a * a + 2.0 * a * (l2 - a) + PI * r * a
            - (under_circle(r + a, a) - under_circle(r - a, a)));
    for (corners, want) in [
        (PipeCorners::Round, round),
        (PipeCorners::Extended, extended),
    ] {
        let mut model = Model::new();
        let spine = arc_then_line(&mut model, r, l2);
        let profile = square_section(
            &mut model,
            Point::new(-r, r, 0.0),
            Vector::new(0.0, -1.0, 0.0),
            a,
        );
        let pipe = make_pipe_shell_with(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            corners,
            1e-3,
            T,
        )
        .unwrap_or_else(|e| panic!("{corners:?}: {e}"))
        .shape;
        let found = check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{corners:?}: {found}");
        let (v, _) = measure(&model, &pipe);
        assert!(
            (v - want).abs() < want * 1e-9,
            "{corners:?}: {v} against {want}"
        );
    }
}

/// The spine: the quarter arc of [`arc_then_line`], then a sharp quarter
/// turn left into a quarter arc of radius `r2` about (-r2, 0), from the
/// origin heading +y to (-r2, r2). Both arcs bend left, so their Frenet
/// frames both have the binormal +z.
fn arc_then_arc(model: &mut Model, r: f64, r2: f64) -> Shape {
    let first = Frame::new(
        Point::new(0.0, r, 0.0),
        Direction::Z,
        Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
        T,
    )
    .unwrap();
    let arc = make_edge(
        model,
        CircleCurve::new(Circle::new(first, r, T).unwrap()).into(),
        (1.5 * PI, 2.0 * PI),
        T,
    )
    .unwrap()
    .shape;
    let (_, joint) = edge_vertices(model, &arc).unwrap().unwrap();
    let second = Frame::new(Point::new(-r2, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let end = make_vertex(model, Point::new(-r2, r2, 0.0)).shape;
    let turn = make_edge_between(
        model,
        CircleCurve::new(Circle::new(second, r2, T).unwrap()).into(),
        (0.0, 0.5 * PI),
        &joint,
        &end,
        T,
    )
    .unwrap()
    .shape;
    make_wire(model, &[arc, turn], T).unwrap().shape
}

/// The area the two arcs' bands share inside the corner: the region over
/// `y > 0` and left of `x = 0` between the two inner circles, the first
/// about (0, r) of radius `r - a` and the second about (-r2, 0) of radius
/// `r2 - a`. Its corners are the origin, (0, a), the circles' crossing `p`
/// and (-a, 0); it is that quadrilateral less the two circular segments the
/// circles cut off it, `ρ² (θ - sin θ) / 2` each.
fn inner_overlap(r: f64, r2: f64, a: f64) -> f64 {
    let (c1, rho1) = (Point::new(0.0, r, 0.0), r - a);
    let (c2, rho2) = (Point::new(-r2, 0.0, 0.0), r2 - a);
    let d = c1.distance(c2);
    let along = (d * d + rho1 * rho1 - rho2 * rho2) / (2.0 * d);
    let off = (rho1 * rho1 - along * along).sqrt();
    let u = (c2 - c1) / d;
    let w = Vector::new(-u.y, u.x, 0.0);
    let p = [c1 + u * along + w * off, c1 + u * along - w * off]
        .into_iter()
        .min_by(|x, y| {
            x.distance(Point::ORIGIN)
                .total_cmp(&y.distance(Point::ORIGIN))
        })
        .unwrap();
    let quad = [
        Point::ORIGIN,
        Point::new(0.0, a, 0.0),
        p,
        Point::new(-a, 0.0, 0.0),
    ];
    let area = (0..4)
        .map(|i| {
            let (s, t) = (quad[i], quad[(i + 1) % 4]);
            s.x * t.y - t.x * s.y
        })
        .sum::<f64>()
        .abs()
        / 2.0;
    let segment = |rho: f64, chord: f64| {
        let theta = 2.0 * (chord / (2.0 * rho)).asin();
        rho * rho * (theta - theta.sin()) / 2.0
    };
    area - segment(rho1, quad[1].distance(p)) - segment(rho2, p.distance(quad[3]))
}

/// A square down [`arc_then_arc`] under the Frenet frame. The arcs' bands
/// are annular sectors of area `π r a` and `π r2 a`, sharing
/// [`inner_overlap`] inside the corner.
///
/// Round: the outer half of the section turned about the corner adds a
/// quarter disc of radius `a` in `x > 0, y < 0`, which neither band
/// reaches.
///
/// Extended: the first arc runs on along +x over `[0, a] × [-a, a]` and the
/// second back along -y over `[-a, a] × [-a, 0]`, the two sharing `a²`.
/// The first run-on overlaps the second band where `x < √((r2 + a)² - y²) -
/// r2`, `∫₀ᵃ √((r2 + a)² - y²) dy - r2 a`, and the second the first band
/// likewise with `r`.
#[test]
fn a_square_under_the_frenet_frame_turns_a_corner_between_arcs_each_way() {
    let (a, r, r2) = (1.0, 8.0, 5.0);
    let bands = PI * r * a + PI * r2 * a - inner_overlap(r, r2, a);
    let round = 2.0 * a * (bands + PI * a * a / 4.0);
    let extended = 2.0
        * a
        * (bands + 4.0 * a * a
            - a * a
            - (under_circle(r2 + a, a) - r2 * a)
            - (under_circle(r + a, a) - r * a));
    for (corners, want) in [
        (PipeCorners::Round, round),
        (PipeCorners::Extended, extended),
    ] {
        let mut model = Model::new();
        let spine = arc_then_arc(&mut model, r, r2);
        let profile = square_section(
            &mut model,
            Point::new(-r, r, 0.0),
            Vector::new(0.0, -1.0, 0.0),
            a,
        );
        let pipe = make_pipe_shell_with(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::Frenet,
            corners,
            1e-3,
            T,
        )
        .unwrap_or_else(|e| panic!("{corners:?}: {e}"))
        .shape;
        let found = check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{corners:?}: {found}");
        let (v, _) = measure(&model, &pipe);
        assert!(
            (v - want).abs() < want * 1e-9,
            "{corners:?}: {v} against {want}"
        );
    }
}

/// The square ring of side `l` under the rotation-minimizing frame, which
/// comes home untwisted round a plane ring. Extended, every corner's
/// square `[-a, a]²` is covered whole and the bands meet nowhere else: the
/// bands' area `8 a l`. Round, each corner covers three quarters of its
/// square and a quarter disc: `a² - π a² / 4` less per corner.
#[test]
fn a_square_round_a_square_ring_turns_its_corners_each_way() {
    let (a, l) = (1.0, 10.0);
    let extended = 2.0 * a * 8.0 * a * l;
    let round = 2.0 * a * (8.0 * a * l - 4.0 * (a * a - PI * a * a / 4.0));
    for (corners, want) in [
        (PipeCorners::Extended, extended),
        (PipeCorners::Round, round),
    ] {
        let mut model = Model::new();
        let profile = square_section(&mut model, Point::ORIGIN, Vector::new(1.0, 0.0, 0.0), a);
        let spine = make_polygon(
            &mut model,
            &[
                Point::ORIGIN,
                Point::new(l, 0.0, 0.0),
                Point::new(l, l, 0.0),
                Point::new(0.0, l, 0.0),
            ],
            true,
            T,
        )
        .unwrap()
        .shape;
        let pipe = make_pipe_shell_with(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            corners,
            1e-3,
            T,
        )
        .unwrap_or_else(|e| panic!("{corners:?}: {e}"))
        .shape;
        let found = check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{corners:?}: {found}");
        let (v, _) = measure(&model, &pipe);
        assert!(
            (v - want).abs() < want * 1e-9,
            "{corners:?}: {v} against {want}"
        );
    }
}

/// A disc round a circle under the binormal along the circle's axis: the
/// rotation-minimizing frame, so the torus, `2 π R π r²` by Pappus. Tilted
/// off the axis, the binormal's frames are not, and the disc is skinned
/// round them: the same torus, held to the skin's tolerance and the mesh's
/// chord over the wall's area.
#[test]
fn a_binormal_round_a_circle_is_the_torus_tilted_or_not() {
    let (big, small) = (10.0, 1.0);
    let mut model = Model::new();
    let ring_frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
    let ring = Circle::new(ring_frame, big, T).unwrap();
    let edge = make_edge(
        &mut model,
        CircleCurve::new(ring).into(),
        (0.0, 2.0 * PI),
        T,
    )
    .unwrap()
    .shape;
    let spine = make_wire(&mut model, &[edge], T).unwrap().shape;
    let at = Point::new(big, 0.0, 0.0);
    let disc_frame = Frame::new(at, Direction::Y, Direction::X, T).unwrap();
    let disc_circle = Circle::new(disc_frame, small, T).unwrap();
    let disc_edge = make_edge(
        &mut model,
        CircleCurve::new(disc_circle).into(),
        (0.0, 2.0 * PI),
        T,
    )
    .unwrap()
    .shape;
    let disc_wire = make_wire(&mut model, &[disc_edge], T).unwrap().shape;
    let disc = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(disc_frame)).into(),
        &[disc_wire],
        T,
    )
    .unwrap()
    .shape;
    let torus = make_pipe_shell_with(
        &mut model,
        &disc,
        &spine,
        &PipeLaw::Binormal(Direction::Z),
        PipeCorners::Mitre,
        1e-3,
        T,
    )
    .unwrap()
    .shape;
    assert!(check(&model, &torus, T).unwrap().is_valid());
    let (v, _) = measure(&model, &torus);
    let want = 2.0 * PI * big * PI * small * small;
    assert!((v - want).abs() < want * 1e-6, "{v} against {want}");

    let skinned = make_pipe_shell_with(
        &mut model,
        &disc,
        &spine,
        &tilted(),
        PipeCorners::Mitre,
        1e-3,
        T,
    )
    .unwrap()
    .shape;
    assert!(check(&model, &skinned, T).unwrap().is_valid());
    let (v, _) = measure(&model, &skinned);
    let bound = (1e-3 + 1e-4) * 2.0 * PI * small * 2.0 * PI * (big + small);
    assert!(
        (v - want).abs() < bound,
        "{v} against {want}, bound {bound}"
    );
}

#[test]
fn corners_the_law_cannot_turn_are_refused_by_name() {
    let (a, l) = (1.0, 10.0);
    let mut model = Model::new();
    // A rectangle: the tilted binormal's section on the second leg is not
    // the first leg's turned about the corner, and no symmetry hides that.
    let profile = rectangle_section(
        &mut model,
        Point::ORIGIN,
        Vector::new(0.0, 1.0, 0.0),
        a,
        2.0 * a,
    );
    let spine = make_polygon(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(0.0, l, 0.0),
            Point::new(l, l, 0.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let err = make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &tilted(),
        PipeCorners::Round,
        1e-3,
        T,
    )
    .unwrap_err();
    assert!(
        err.to_string().contains("a round corner turns the section"),
        "{err}"
    );
    let err = make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::Frenet,
        PipeCorners::Extended,
        1e-3,
        T,
    )
    .unwrap_err();
    assert!(
        err.to_string()
            .contains("the Frenet frame is not defined along a straight run"),
        "{err}"
    );
}
