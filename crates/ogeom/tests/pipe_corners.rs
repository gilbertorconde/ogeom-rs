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

/// A circle of `radius` about `at`, square to `heading` (in the `xy`
/// plane), as a wire.
fn circle_ring(model: &mut Model, at: Point, heading: Vector, radius: f64) -> Shape {
    let across = Vector::new(0.0, 0.0, 1.0).cross(heading);
    let frame = Frame::new(
        at,
        Direction::new(heading, T).unwrap(),
        Direction::new(across, T).unwrap(),
        T,
    )
    .unwrap();
    let circle = Circle::new(frame, radius, T).unwrap();
    let edge = make_edge(model, CircleCurve::new(circle).into(), (0.0, 2.0 * PI), T)
        .unwrap()
        .shape;
    make_wire(model, &[edge], T).unwrap().shape
}

/// A face on the plane through the origin square to `y`, bounded by
/// `rings`, the outer first.
fn face_of(model: &mut Model, rings: &[Shape]) -> Shape {
    let frame = Frame::new(Point::ORIGIN, Direction::Y, Direction::X, T).unwrap();
    make_face(model, PlaneSurface::new(Plane::new(frame)).into(), rings, T)
        .unwrap()
        .shape
}

/// An auxiliary spine at height `h` beside the L (0,0,0) -> (0,L,0) ->
/// (L,L,0), its legs off the L's so the section's normal, pointing at it,
/// turns along both of the L's legs. Each leg's planes cross it once:
/// the first leg's on the part left of the L, the second's on the part
/// beyond it. Where `still`, it runs parallel to each of the L's legs for
/// the last 5 before the corner and the first 5 after it, so the frame
/// holds still there, but each leg's frame differs from the other's
/// mirrored across the corner.
fn guide_beside_an_l(model: &mut Model, l: f64, h: f64, still: bool) -> Shape {
    let points: Vec<Point> = if still {
        vec![
            Point::new(-3.0, -1.0, h),
            Point::new(-1.0, l - 5.0, h),
            Point::new(-1.0, l + 2.0, h),
            Point::new(5.0, l + 2.0, h),
            Point::new(l + 1.0, l + 5.0, h),
        ]
    } else {
        vec![
            Point::new(-3.0, -1.0, h),
            Point::new(-1.0, l + 1.0, h),
            Point::new(l + 1.0, l + 4.0, h),
        ]
    };
    make_polygon(model, &points, false, T).unwrap().shape
}

/// A rectangle `2 a` by `4 a` with a circular hole of radius `r` about
/// its centre, and a disc of radius `2 a` with that hole, down the L
/// under an auxiliary spine that turns the section along both legs, so
/// each leg is a skinned run and the two meet at the mitre with their
/// frames there unrelated. The hole's circle is centred on the spine and
/// square to it, so every frame carries it onto itself: its wall down
/// each leg is the cylinder, cut at the mitre in one ellipse.
///
/// Where the frame holds still within the section's reach of the corner,
/// each run there is the prism of its section, and the mitred L holds the
/// section's area times the spine's length (the wedge each prism runs on
/// past the corner and the part the mitre plane cuts back balance about
/// the section's centre): within the skin's tolerance over the outer
/// ring's walls. Where the frame turns up to the corner, the outer ring's
/// turned prisms no longer balance, but the hole's tube still does, and
/// the holed pipe is the rectangle's under the same law less `π r² 2 L`.
/// The annulus's two rings are both such circles, so it is the mitred
/// tube under any frame, `π (4 a² - r²) 2 L`.
#[test]
fn a_holed_section_down_an_l_mitres_under_a_turning_auxiliary_spine() {
    let (a, l, r, h, skin) = (1.0, 10.0, 0.5, 5.0, 1e-3);
    let y = Vector::new(0.0, 1.0, 0.0);
    let hole_volume = PI * r * r * 2.0 * l;
    let sweep = |section: &str, still: bool| {
        let mut model = Model::new();
        let rect = rectangle_section(&mut model, Point::ORIGIN, y, a, 2.0 * a);
        let rect = model.ordered_children_of(&rect).unwrap()[0].clone();
        let hole = circle_ring(&mut model, Point::ORIGIN, y, r);
        let profile = match section {
            "holed" => face_of(&mut model, &[rect, hole]),
            "solid" => face_of(&mut model, &[rect]),
            _ => {
                let disc = circle_ring(&mut model, Point::ORIGIN, y, 2.0 * a);
                face_of(&mut model, &[disc, hole])
            }
        };
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
        let guide = guide_beside_an_l(&mut model, l, h, still);
        let pipe = make_pipe_shell_with(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::Auxiliary { guide: &guide },
            PipeCorners::Mitre,
            skin,
            T,
        )
        .unwrap_or_else(|e| panic!("{section}, still {still}: {e}"))
        .shape;
        let found = check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{section}, still {still}: {found}");
        let (v, exact) = measure(&model, &pipe);
        assert!(exact, "{section}, still {still}: not integrated exactly");
        v
    };

    let v = sweep("holed", true);
    let want = 8.0 * a * a * 2.0 * l - hole_volume;
    let reach = (a * a + 4.0 * a * a).sqrt();
    let bound = skin * 12.0 * a * (2.0 * l + 2.0 * reach);
    assert!(
        (v - want).abs() < bound,
        "{v} against {want}, bound {bound}"
    );

    let (holed, solid) = (sweep("holed", false), sweep("solid", false));
    assert!(
        (solid - holed - hole_volume).abs() < hole_volume * 1e-9,
        "{solid} less {holed} against the tube's {hole_volume}"
    );

    let v = sweep("annulus", false);
    let want = PI * (4.0 * a * a - r * r) * 2.0 * l;
    assert!((v - want).abs() < want * 1e-12, "{v} against {want}");
}

/// An auxiliary spine at height `h` beside the L (0,0,0) -> (0,L,0) ->
/// (L,L,0) that each station's plane crosses on a line parallel to the
/// leg there: at `x = 2` for the first 2 of the first leg, at `x = -1`
/// for its last 3, and at `y = L + 2` for the first 3 of the second leg.
/// The frame holds still on each of those stretches, its normal there
/// along `(2, 0, h)`, `(-1, 0, h)` and `(0, 2, h)`, and turns between, so
/// both legs are skinned runs.
fn guide_held_at_the_ends(model: &mut Model, l: f64, h: f64) -> Shape {
    let points = [
        Point::new(2.0, -1.0, h),
        Point::new(2.0, 2.0, h),
        Point::new(-1.0, l - 3.0, h),
        Point::new(-1.0, l + 2.0, h),
        Point::new(3.0, l + 2.0, h),
        Point::new(l + 1.0, l + 5.0, h),
    ];
    make_polygon(model, &points, false, T).unwrap().shape
}

/// Sections off the spine down the L under [`guide_held_at_the_ends`]:
/// a disc, a rectangle, and a rectangle about the spine with a hole off
/// it, each meeting the mitre between two skinned runs.
///
/// A turn about a straight spine moves no volume, so each run holds its
/// section's area `A` times its length wherever the frame turns. Within
/// the section's reach of the corner the frame holds still and each run
/// is the prism of its section, run on to the mitre plane or cut back by
/// it: a point `q` off the corner runs `-q . d₂` past the first leg's end
/// and `q . d₁` back from the second leg's start, `d₁` and `d₂` the legs'
/// headings. So the pipe holds `A (2 L + c₂ . d₁ - c₁ . d₂)`, `c₁` and
/// `c₂` the section's centroid off the corner in each leg's frame there:
/// the start frame's normal and binormal components carried to each.
/// A hole subtracts its own such term. The volume strays from that by at
/// most the skin's tolerance over the walls: the perimeter times the
/// longest path a point runs, the spine and a reach past the corner.
#[test]
fn an_off_centre_section_down_an_l_mitres_between_skinned_runs() {
    let (l, h, skin) = (10.0, 5.0, 1e-3);
    let unit = |v: Vector| v / v.magnitude();
    let (d1, d2) = (Vector::new(0.0, 1.0, 0.0), Vector::new(1.0, 0.0, 0.0));
    let across = Vector::new(0.0, 0.0, 1.0).cross(d1);
    let up = Vector::new(0.0, 0.0, 1.0);
    let n0 = unit(Vector::new(2.0, 0.0, h));
    let n1 = unit(Vector::new(-1.0, 0.0, h));
    let n2 = unit(Vector::new(0.0, 2.0, h));
    // The mitre's share for a section whose centroid stands `s` across and
    // `t` up from the spine at its start.
    let mitred = |s: f64, t: f64| {
        let c = across * s + up * t;
        let (u, v) = (c.dot(n0), c.dot(d1.cross(n0)));
        let c1 = n1 * u + d1.cross(n1) * v;
        let c2 = n2 * u + d2.cross(n2) * v;
        c2.dot(d1) - c1.dot(d2)
    };
    let (s, t, a, r) = (0.6, -0.8, 1.0, 0.5);
    let at = |s: f64, t: f64| Point::ORIGIN + across * s + up * t;
    let rect = |model: &mut Model, s: f64, t: f64| {
        let pts = [(-a, -2.0 * a), (a, -2.0 * a), (a, 2.0 * a), (-a, 2.0 * a)]
            .map(|(x, y)| at(s + x, t + y));
        make_polygon(model, &pts, true, T).unwrap().shape
    };
    let sweep = |section: &str| {
        let mut model = Model::new();
        let rings = match section {
            "disc" => vec![circle_ring(&mut model, at(s, t), d1, a)],
            "rectangle" => vec![rect(&mut model, s, t)],
            _ => {
                let outer = rect(&mut model, 0.0, 0.0);
                let hole = circle_ring(&mut model, at(s / 2.0, t), d1, r);
                vec![outer, hole]
            }
        };
        let profile = face_of(&mut model, &rings);
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
        let guide = guide_held_at_the_ends(&mut model, l, h);
        let pipe = make_pipe_shell_with(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::Auxiliary { guide: &guide },
            PipeCorners::Mitre,
            skin,
            T,
        )
        .unwrap_or_else(|e| panic!("{section}: {e}"))
        .shape;
        let found = check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{section}: {found}");
        let (v, exact) = measure(&model, &pipe);
        assert!(exact, "{section}: not integrated exactly");
        v
    };
    let path = 2.0 * l + 2.0 * (s.hypot(t) + (5.0 * a * a).sqrt());
    for (section, want, perimeter) in [
        ("disc", PI * a * a * (2.0 * l + mitred(s, t)), 2.0 * PI * a),
        (
            "rectangle",
            8.0 * a * a * (2.0 * l + mitred(s, t)),
            12.0 * a,
        ),
        (
            "holed",
            8.0 * a * a * 2.0 * l - PI * r * r * (2.0 * l + mitred(s / 2.0, t)),
            12.0 * a + 2.0 * PI * r,
        ),
    ] {
        let v = sweep(section);
        let bound = skin * perimeter * path;
        assert!(
            (v - want).abs() < bound,
            "{section}: {v} against {want}, bound {bound}"
        );
    }
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
