//! Open wires swept as sheets of fitted skins: down spines with a
//! free-form leg, and under laws turning the section other than the
//! rotation-minimizing frame does, against closed forms and against the
//! walls of the face the open wires halve.
//!
//! Down a planar spine the rotation-minimizing frame keeps the section's
//! frame in the spine's plane: a profile point at offset `x` along the
//! in-plane normal `e` (square to the heading, `+x` where the spine leaves
//! the origin along `+z`) runs a path of length `L - x Θ` down a curved
//! leg of length `L` whose heading turns through `Θ` toward `e`, wherever
//! `x` times the curvature stays under one. A corner turning through `θ`
//! toward `+x` has its outside at `α = -x`: mitred, a point outside gains
//! `2 α tan(θ/2)` on its path; rounded, `α θ`. A sheet's area is its
//! paths' lengths integrated along the profile; its boundary is the
//! profile twice and the paths of its two ends.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_algo::{check, curve_length, make_polygon, surface_properties};
use ogeom_core::{FastMap, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, KnotVector, Point, Vector};
use ogeom_mesh::Deflection;
use ogeom_offset::{PipeCorners, PipeLaw, make_pipe_shell_with};
use ogeom_topo::{EdgeRepr, Filter, Model, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

/// The skins' tolerance.
const SKIN: f64 = 1e-5;

/// Whether a measure of a sheet whose walls stand within `SKIN` of the
/// exact ones meets its closed form: moving a wall by `ε` moves its area by
/// at most `ε` times its boundary's length and its area times its greatest
/// curvature, under one on these sheets, and its boundary's length by
/// `ε` along each of its curves, a few dozen at most.
fn near(got: f64, want: f64, boundary: f64, area: f64) -> bool {
    (got - want).abs() <= SKIN * (boundary + area)
}

fn sweep_with(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    corners: PipeCorners,
) -> Shape {
    make_pipe_shell_with(model, profile, spine, law, corners, SKIN, T)
        .unwrap_or_else(|e| panic!("{corners:?}: {e}"))
        .shape
}

fn sweep(model: &mut Model, profile: &Shape, spine: &Shape, corners: PipeCorners) -> Shape {
    sweep_with(model, profile, spine, &PipeLaw::RotationMinimizing, corners)
}

fn refusal(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    corners: PipeCorners,
) -> String {
    match make_pipe_shell_with(model, profile, spine, law, corners, SKIN, T) {
        Ok(_) => panic!("{corners:?}: built"),
        Err(e) => e.to_string(),
    }
}

/// A sheet's area and the length of its boundary, the edges of one face,
/// after checking it is one shell with no finding but its being open.
fn measured(model: &Model, sheet: &Shape) -> (f64, f64) {
    assert_eq!(model.kind_of(sheet).unwrap(), ShapeType::Shell);
    let diagnosis = check(model, sheet, T).unwrap();
    assert!(
        diagnosis
            .problems
            .iter()
            .all(|p| p.kind == ShapeType::Shell),
        "{diagnosis}"
    );
    let area = surface_properties(model, sheet, Deflection::with_chord(1e-4).unwrap(), T)
        .unwrap()
        .mass;
    let mut uses: FastMap<ogeom_topo::TShapeId, (Shape, usize)> = FastMap::default();
    for face in explore(model, sheet, Filter::OfType(ShapeType::Face)).unwrap() {
        for edge in explore(model, &face, Filter::OfType(ShapeType::Edge)).unwrap() {
            uses.entry(edge.node()).or_insert((edge, 0)).1 += 1;
        }
    }
    assert!(
        uses.values().all(|(_, n)| *n <= 2),
        "an edge of three faces"
    );
    let mut boundary = 0.0;
    for (edge, _) in uses.values().filter(|(_, n)| *n == 1) {
        let data = model.node(edge).unwrap().data().as_edge().unwrap();
        if data.degenerate {
            continue;
        }
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            panic!("an edge with no curve");
        };
        let curve = model.geometry().curve(*curve).unwrap();
        boundary += curve_length(curve, *range, T).unwrap();
    }
    (area, boundary)
}

/// An open polyline in the plane `z = 0`, square to a spine leaving the
/// origin along `+z`.
fn polyline(model: &mut Model, points: &[(f64, f64)]) -> Shape {
    let points: Vec<Point> = points.iter().map(|&(x, y)| Point::new(x, y, 0.0)).collect();
    make_polygon(model, &points, false, T).unwrap().shape
}

/// A cubic Bézier in the plane `y = 0` from `a` to `b`, its inner control
/// points `p` and `q`.
fn bezier(a: Point, p: Point, q: Point, b: Point) -> ogeom_geom::Curve {
    ogeom_geom::Curve::BSpline(
        ogeom_geom::BSplineCurve::new(
            KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap(),
            vec![a, p, q, b],
            T,
        )
        .unwrap(),
    )
}

/// A wire through `curves`, each starting where the last ends.
fn chained(model: &mut Model, curves: &[ogeom_geom::Curve]) -> Shape {
    let mut edges = Vec::new();
    let first = curves[0].point_at(curves[0].domain().0, T).unwrap();
    let mut from = ogeom_algo::make_vertex(model, first).shape;
    for curve in curves {
        let range = curve.domain();
        let to = ogeom_algo::make_vertex(model, curve.point_at(range.1, T).unwrap()).shape;
        edges.push(
            ogeom_algo::make_edge_between(model, curve.clone(), range, &from, &to, T)
                .unwrap()
                .shape,
        );
        from = to;
    }
    ogeom_algo::make_wire(model, &edges, T).unwrap().shape
}

fn line(a: Point, b: Point) -> ogeom_geom::Curve {
    ogeom_geom::LineCurve::segment(a, b, T).unwrap().into()
}

/// The length of a curve and the angle its heading turns through about
/// `+y`, from `+z` toward `+x`.
fn length_and_turn(curve: &ogeom_geom::Curve) -> (f64, f64) {
    let range = curve.domain();
    let (t0, t1) = (
        curve.d1_at(range.0, T).unwrap(),
        curve.d1_at(range.1, T).unwrap(),
    );
    let turn = t0
        .cross(t1)
        .dot(Vector::new(0.0, 1.0, 0.0))
        .atan2(t0.dot(t1));
    (curve_length(curve, range, T).unwrap(), turn)
}

/// The area and boundary of a polyline profile's sheet where the path of
/// the point `(x, y)` has length `run(x)`, linear in `x`.
fn linear_closed_form(profile: &[(f64, f64)], run: impl Fn(f64) -> f64) -> (f64, f64) {
    let mut length = 0.0;
    let mut area = 0.0;
    for pair in profile.windows(2) {
        let s = ((pair[1].0 - pair[0].0).powi(2) + (pair[1].1 - pair[0].1).powi(2)).sqrt();
        length += s;
        area += s * run(f64::midpoint(pair[0].0, pair[1].0));
    }
    let rails = run(profile[0].0) + run(profile[profile.len() - 1].0);
    (area, 2.0 * length + rails)
}

/// A free-form leg bending toward `+x`, then back: an S whose heading
/// turns through a net angle, in the plane `y = 0`.
fn s_leg(from: Point, heading: Vector) -> ogeom_geom::Curve {
    let side = Vector::new(heading.z, 0.0, -heading.x);
    bezier(
        from,
        from + heading * 4.0,
        from + heading * 6.0 + side * 3.0,
        from + heading * 9.0 + side * 3.5,
    )
}

#[test]
fn an_open_profile_down_a_planar_spline_sweeps_to_its_closed_form() {
    let leg = s_leg(Point::ORIGIN, Vector::new(0.0, 0.0, 1.0));
    let (length, turn) = length_and_turn(&leg);
    for profile in [
        vec![(-1.0, 0.5), (1.5, 0.5)],
        vec![(0.5, -1.0), (0.5, 1.0)],
        vec![(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0)],
    ] {
        for corners in [PipeCorners::Mitre, PipeCorners::Round] {
            let mut model = Model::new();
            let spine = chained(&mut model, std::slice::from_ref(&leg));
            let wire = polyline(&mut model, &profile);
            let sheet = sweep(&mut model, &wire, &spine, corners);
            let (area, boundary) = measured(&model, &sheet);
            let (want_area, want_boundary) = linear_closed_form(&profile, |x| length - x * turn);
            assert!(
                near(area, want_area, want_boundary, want_area)
                    && near(boundary, want_boundary, want_boundary, want_area),
                "{profile:?} {corners:?}: area {area} against {want_area}, \
                 boundary {boundary} against {want_boundary}"
            );
        }
    }
}

/// A half circle of radius one about the spine's start, from `(1, 0)`
/// round through `(0, 1)` to `(-1, 0)` where `turned`, else from `(0, -1)`
/// through `(1, 0)` to `(0, 1)`.
fn half_circle(model: &mut Model, turned: bool) -> Shape {
    let frame = ogeom_math::Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
    let circle = ogeom_math::Circle::new(frame, 1.0, T).unwrap();
    let range = if turned {
        (0.0, PI)
    } else {
        (-PI / 2.0, PI / 2.0)
    };
    let edge = ogeom_algo::make_edge(model, ogeom_geom::CircleCurve::new(circle).into(), range, T)
        .unwrap()
        .shape;
    ogeom_algo::make_wire(model, &[edge], T).unwrap().shape
}

/// An arc of the profile sweeps its length times the leg less `Θ` times
/// the integral of `x` along it: two for the half circle on `x ≥ 0`.
#[test]
fn a_curved_profile_edge_down_a_planar_spline_sweeps_to_its_closed_form() {
    let leg = s_leg(Point::ORIGIN, Vector::new(0.0, 0.0, 1.0));
    let (length, turn) = length_and_turn(&leg);
    let mut model = Model::new();
    let spine = chained(&mut model, std::slice::from_ref(&leg));
    let wire = half_circle(&mut model, false);
    let sheet = sweep(&mut model, &wire, &spine, PipeCorners::Mitre);
    let (area, boundary) = measured(&model, &sheet);
    let (want_area, want_boundary) = (PI * length - 2.0 * turn, 2.0 * PI + 2.0 * length);
    assert!(
        near(area, want_area, want_boundary, want_area)
            && near(boundary, want_boundary, want_boundary, want_area),
        "area {area} against {want_area}, boundary {boundary} against {want_boundary}"
    );
}

#[test]
fn an_open_profile_outside_a_corner_against_a_spline_sweeps_to_its_closed_form() {
    let theta = PI / 3.0;
    let corner = Point::new(0.0, 0.0, 5.0);
    let heading = Vector::new(theta.sin(), 0.0, theta.cos());
    let leg = s_leg(corner, heading);
    let (length, turn) = length_and_turn(&leg);
    let profile = [(-3.0, -1.0), (-1.0, -1.0), (-1.0, 1.0), (0.0, 1.0)];
    for (curves, before) in [
        (vec![line(Point::ORIGIN, corner), leg.clone()], 5.0),
        // A spline into the corner: straight, so its own closed form.
        (
            vec![
                bezier(
                    Point::ORIGIN,
                    Point::new(0.0, 0.0, 1.0),
                    Point::new(0.0, 0.0, 3.0),
                    corner,
                ),
                leg.clone(),
            ],
            5.0,
        ),
    ] {
        for corners in [PipeCorners::Mitre, PipeCorners::Round] {
            let mut model = Model::new();
            let spine = chained(&mut model, &curves);
            let wire = polyline(&mut model, &profile);
            let sheet = sweep(&mut model, &wire, &spine, corners);
            let (area, boundary) = measured(&model, &sheet);
            let gain = |x: f64| {
                let alpha = -x;
                if corners == PipeCorners::Mitre {
                    2.0 * alpha * (theta / 2.0).tan()
                } else {
                    alpha * theta
                }
            };
            let (want_area, want_boundary) =
                linear_closed_form(&profile, |x| before + length - x * turn + gain(x));
            assert!(
                near(area, want_area, want_boundary, want_area)
                    && near(boundary, want_boundary, want_boundary, want_area),
                "{corners:?}: area {area} against {want_area}, \
                 boundary {boundary} against {want_boundary}"
            );
        }
    }
}

/// The face swept with the profile's closed wire bounded: its area less
/// its two caps, `cap` each.
fn face_walls(
    model: &mut Model,
    closed: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    cap: f64,
    skin: f64,
) -> f64 {
    let face = ogeom_algo::make_face(
        model,
        ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(
            ogeom_math::Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap(),
        ))
        .into(),
        std::slice::from_ref(closed),
        T,
    )
    .unwrap()
    .shape;
    let solid = make_pipe_shell_with(model, &face, spine, law, PipeCorners::Mitre, skin, T)
        .unwrap_or_else(|e| panic!("the face: {e}"))
        .shape;
    surface_properties(model, &solid, Deflection::with_chord(1e-4).unwrap(), T)
        .unwrap()
        .mass
        - 2.0 * cap
}

/// Inside a corner against a curved leg, a mitred face ends each leg's
/// walls where its points' paths cross, as the sheet does: the two open
/// halves of a rectangle across the corner's axis sum to the face's walls,
/// within what both skins are held to. Rounded, the halves gain on the
/// mitred ones only outside the turn, `∫ α (θ - 2 tan(θ/2))` along the
/// profile.
#[test]
fn two_open_halves_across_a_corner_against_a_spline_bound_its_face() {
    let theta = PI / 3.0;
    let corner = Point::new(0.0, 0.0, 5.0);
    let heading = Vector::new(theta.sin(), 0.0, theta.cos());
    let leg = s_leg(corner, heading);
    let rectangle = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
    let halves = [
        vec![rectangle[0], rectangle[1], rectangle[2]],
        vec![rectangle[2], rectangle[3], rectangle[0]],
    ];
    let mut model = Model::new();
    let spine = chained(&mut model, &[line(Point::ORIGIN, corner), leg.clone()]);
    let mut sums = [0.0; 2];
    let mut boundary = 0.0;
    for (sum, corners) in sums
        .iter_mut()
        .zip([PipeCorners::Mitre, PipeCorners::Round])
    {
        for half in &halves {
            let wire = polyline(&mut model, half);
            let sheet = sweep(&mut model, &wire, &spine, corners);
            let (area, edges) = measured(&model, &sheet);
            *sum += area;
            boundary = edges;
        }
    }
    let points: Vec<Point> = rectangle
        .iter()
        .map(|&(x, y)| Point::new(x, y, 0.0))
        .collect();
    let closed = make_polygon(&mut model, &points, true, T).unwrap().shape;
    // The face's skin is held to a tenth of a micrometre.
    let face_skin = 1e-4;
    let walls = face_walls(
        &mut model,
        &closed,
        &spine,
        &PipeLaw::RotationMinimizing,
        4.0,
        face_skin,
    );
    let bound = (SKIN + face_skin) * (2.0 * boundary + walls);
    assert!(
        (sums[0] - walls).abs() <= bound,
        "mitred halves {} against the face's walls {walls}",
        sums[0]
    );
    // `α = -x` is one along the side at `x = -1` and runs to naught along
    // the halves of the two sides across the axis: three in all.
    let outside = 3.0;
    let want = outside * (theta - 2.0 * (theta / 2.0).tan());
    assert!(
        (sums[1] - sums[0] - want).abs() <= SKIN * (2.0 * boundary + walls),
        "rounded halves gain {} on mitred ones against {want}",
        sums[1] - sums[0]
    );
}

/// A curved profile edge across a rounded corner's axis is cut there; its
/// sheet gains on the mitred one only outside the turn, `∫ α` being one
/// along the quarter of the half circle there. Mitred, the half circle and
/// its diameter sum to the walls of the half disc they bound.
#[test]
fn a_curved_profile_edge_across_a_corner_against_a_spline_is_cut_on_its_axis() {
    let theta = PI / 3.0;
    let corner = Point::new(0.0, 0.0, 5.0);
    let heading = Vector::new(theta.sin(), 0.0, theta.cos());
    let leg = s_leg(corner, heading);
    let mut model = Model::new();
    let spine = chained(&mut model, &[line(Point::ORIGIN, corner), leg.clone()]);
    let arc = half_circle(&mut model, true);
    let mut areas = [0.0; 2];
    let mut boundary = 0.0;
    for (area, corners) in areas
        .iter_mut()
        .zip([PipeCorners::Mitre, PipeCorners::Round])
    {
        let sheet = sweep(&mut model, &arc, &spine, corners);
        (*area, boundary) = measured(&model, &sheet);
    }
    let want = theta - 2.0 * (theta / 2.0).tan();
    assert!(
        (areas[1] - areas[0] - want).abs() <= SKIN * (2.0 * boundary + areas[0]),
        "the rounded half circle gains {} on the mitred one against {want}",
        areas[1] - areas[0]
    );
    let diameter = polyline(&mut model, &[(-1.0, 0.0), (1.0, 0.0)]);
    let sheet = sweep(&mut model, &diameter, &spine, PipeCorners::Mitre);
    let (straight, _) = measured(&model, &sheet);
    let closed = {
        let (a, b) = (Point::new(1.0, 0.0, 0.0), Point::new(-1.0, 0.0, 0.0));
        let frame = ogeom_math::Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
        let circle: ogeom_geom::Curve =
            ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(frame, 1.0, T).unwrap()).into();
        let (va, vb) = (
            ogeom_algo::make_vertex(&mut model, a).shape,
            ogeom_algo::make_vertex(&mut model, b).shape,
        );
        let edges = [
            ogeom_algo::make_edge_between(&mut model, circle, (0.0, PI), &va, &vb, T)
                .unwrap()
                .shape,
            ogeom_algo::make_edge_between(&mut model, line(b, a), (0.0, 2.0), &vb, &va, T)
                .unwrap()
                .shape,
        ];
        ogeom_algo::make_wire(&mut model, &edges, T).unwrap().shape
    };
    let face_skin = 1e-4;
    let walls = face_walls(
        &mut model,
        &closed,
        &spine,
        &PipeLaw::RotationMinimizing,
        PI / 2.0,
        face_skin,
    );
    let bound = (SKIN + face_skin) * (2.0 * boundary + walls);
    assert!(
        (areas[0] + straight - walls).abs() <= bound,
        "the mitred half circle and diameter {} against the face's walls {walls}",
        areas[0] + straight
    );
}

/// A spine of a line and the arc it runs on into tangent, turning toward
/// `+x` in the plane `y = 0`.
fn line_and_arc(model: &mut Model) -> Shape {
    let frame = ogeom_math::Frame::new(
        Point::new(4.0, 0.0, 5.0),
        Direction::Y,
        Direction::new(Vector::new(-1.0, 0.0, 0.0), T).unwrap(),
        T,
    )
    .unwrap();
    let arc: ogeom_geom::Curve =
        ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(frame, 4.0, T).unwrap()).into();
    let arc = ogeom_geom::TrimmedCurve::new(arc, 0.0, PI / 2.0, T)
        .unwrap()
        .into();
    chained(
        model,
        &[line(Point::ORIGIN, Point::new(0.0, 0.0, 5.0)), arc],
    )
}

/// A law turning the section other than the rotation-minimizing frame
/// does sweeps each half of a square as the face sweeps the square's
/// walls, under the same law, down a smooth spine of lines and arcs and
/// down a free-form one; the halves sum to those walls within what both
/// skins are held to.
#[test]
fn two_open_halves_under_a_leaning_binormal_bound_the_face_it_sweeps() {
    let lean = Direction::new(Vector::new(0.4, 1.0, 0.0), T).unwrap();
    let square = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
    let halves = [
        vec![square[0], square[1], square[2]],
        vec![square[2], square[3], square[0]],
    ];
    let spines: [fn(&mut Model) -> Shape; 2] = [line_and_arc, |model| {
        let leg = s_leg(Point::ORIGIN, Vector::new(0.0, 0.0, 1.0));
        chained(model, &[leg])
    }];
    for spine in spines {
        let mut model = Model::new();
        let path = spine(&mut model);
        let law = PipeLaw::Binormal(lean);
        let mut sum = 0.0;
        let mut boundary = 0.0;
        for half in &halves {
            let wire = polyline(&mut model, half);
            let sheet = sweep_with(&mut model, &wire, &path, &law, PipeCorners::Mitre);
            let (area, edges) = measured(&model, &sheet);
            sum += area;
            boundary = edges;
        }
        let points: Vec<Point> = square.iter().map(|&(x, y)| Point::new(x, y, 0.0)).collect();
        let closed = make_polygon(&mut model, &points, true, T).unwrap().shape;
        let walls = face_walls(&mut model, &closed, &path, &law, 4.0, SKIN);
        assert!(
            (sum - walls).abs() <= 2.0 * SKIN * (2.0 * boundary + walls),
            "halves {sum} against the face's walls {walls}"
        );
    }
}

/// Round a sharp corner a law whose frames are the rotation-minimizing
/// ones turned by one fixed angle (a binormal square to a planar spine)
/// sweeps the rotation-minimizing frame's sheet; one that leans, or the
/// Frenet frame at a mitre, is refused by name. Down a planar spline that
/// bends one way, the Frenet frame is the rotation-minimizing one turned
/// by a fixed angle too, and sweeps its sheet, its frames carried between
/// stations alike.
#[test]
fn a_law_round_a_corner_against_a_spline_sweeps_as_the_frame_it_turns() {
    let theta = PI / 3.0;
    let corner = Point::new(0.0, 0.0, 5.0);
    let heading = Vector::new(theta.sin(), 0.0, theta.cos());
    let profile = [(-3.0, -1.0), (-1.0, -1.0), (-1.0, 1.0)];
    for corners in [PipeCorners::Mitre, PipeCorners::Round] {
        let mut model = Model::new();
        let spine = chained(
            &mut model,
            &[line(Point::ORIGIN, corner), s_leg(corner, heading)],
        );
        let wire = polyline(&mut model, &profile);
        let carried = sweep(&mut model, &wire, &spine, corners);
        let (want, _) = measured(&model, &carried);
        let law = PipeLaw::Binormal(Direction::Y);
        let sheet = sweep_with(&mut model, &wire, &spine, &law, corners);
        let (area, _) = measured(&model, &sheet);
        assert!(
            (area - want).abs() <= want * 1e-12,
            "{law:?} {corners:?}: {area} against the carried sheet's {want}"
        );
        let lean = PipeLaw::Binormal(Direction::new(Vector::new(0.4, 1.0, 0.0), T).unwrap());
        let message = refusal(&mut model, &wire, &spine, &lean, corners);
        assert!(
            message.contains("the law turns an open wire's section other than"),
            "{corners:?}: {message}"
        );
    }
    let mut model = Model::new();
    let spine = chained(
        &mut model,
        &[line(Point::ORIGIN, corner), s_leg(corner, heading)],
    );
    let wire = polyline(&mut model, &profile);
    let message = refusal(
        &mut model,
        &wire,
        &spine,
        &PipeLaw::Frenet,
        PipeCorners::Mitre,
    );
    assert!(
        message.contains("a Frenet frame has no direction at a corner"),
        "{message}"
    );

    let bend = bezier(
        Point::ORIGIN,
        Point::new(0.0, 0.0, 4.0),
        Point::new(2.0, 0.0, 7.0),
        Point::new(6.0, 0.0, 8.0),
    );
    let mut model = Model::new();
    let spine = chained(&mut model, &[bend]);
    let wire = polyline(&mut model, &profile);
    let carried = sweep(&mut model, &wire, &spine, PipeCorners::Mitre);
    let (want, edges) = measured(&model, &carried);
    let sheet = sweep_with(
        &mut model,
        &wire,
        &spine,
        &PipeLaw::Frenet,
        PipeCorners::Mitre,
    );
    let (area, _) = measured(&model, &sheet);
    assert!(
        near(area, want, 2.0 * edges, 2.0 * want),
        "Frenet {area} against the carried sheet's {want}"
    );
}

/// What no sheet can join is refused by name: a point inside a corner
/// that turns a curved leg out of its plane, whose paths down the two legs
/// pass each other; an extended corner, which squares a solid's outside.
/// The same corner with the profile wholly outside it sweeps.
#[test]
fn open_wire_sheets_down_skinned_legs_it_cannot_join_are_refused_by_name() {
    let corner = Point::new(0.0, 0.0, 5.0);
    // Turned toward `+y`, then bending toward `+x`: out of the turn's plane.
    let heading = Vector::new(0.0, (PI / 3.0).sin(), (PI / 3.0).cos());
    let leg = bezier(
        corner,
        corner + heading * 4.0,
        corner + heading * 6.0 + Vector::new(3.0, 0.0, 0.0),
        corner + heading * 9.0 + Vector::new(3.5, 0.0, 0.0),
    );
    let mut model = Model::new();
    let spine = chained(&mut model, &[line(Point::ORIGIN, corner), leg]);
    let law = PipeLaw::RotationMinimizing;
    for corners in [PipeCorners::Mitre, PipeCorners::Round] {
        let inside = polyline(&mut model, &[(-1.0, 0.5), (1.0, 0.5)]);
        let message = refusal(&mut model, &inside, &spine, &law, corners);
        assert!(
            message.contains("pass each other") && message.contains("without meeting"),
            "{corners:?}: {message}"
        );
        let outside = polyline(&mut model, &[(-1.0, -0.5), (1.0, -0.5)]);
        let sheet = sweep(&mut model, &outside, &spine, corners);
        measured(&model, &sheet);
    }
    let wire = polyline(&mut model, &[(-1.0, -0.5), (1.0, -0.5)]);
    let message = refusal(&mut model, &wire, &spine, &law, PipeCorners::Extended);
    assert!(
        message.contains("an extended corner squares the outside of a solid's section"),
        "{message}"
    );
}

/// Every face of the sheet faces as the prism of the profile along the
/// spine's start heading does, forward and reversed, round mitred and
/// rounded corners: its normal along the profile's heading crossed with
/// the spine's, the profile's segment square to the spine's plane.
#[test]
fn an_open_sheet_down_skinned_legs_faces_as_its_prism() {
    let theta = PI / 3.0;
    let corner = Point::new(0.0, 0.0, 5.0);
    let heading = Vector::new(theta.sin(), 0.0, theta.cos());
    let leg = s_leg(corner, heading);
    for reversed in [false, true] {
        for corners in [PipeCorners::Mitre, PipeCorners::Round] {
            let mut model = Model::new();
            let spine = chained(&mut model, &[line(Point::ORIGIN, corner), leg.clone()]);
            let mut wire = polyline(&mut model, &[(-1.0, -1.0), (-1.0, 1.0)]);
            if reversed {
                wire = wire.reversed();
            }
            let along = if reversed { -1.0 } else { 1.0 };
            let sheet = sweep(&mut model, &wire, &spine, corners);
            let faces = explore(&model, &sheet, Filter::OfType(ShapeType::Face)).unwrap();
            assert!(faces.len() >= 3, "{reversed} {corners:?}");
            for face in &faces {
                let (at, n) = ogeom_algo::face_normal(&model, face, T).unwrap();
                // The sheet's normal lies in the spine's plane, square to
                // the profile's heading `±y`; `+y` crossed with the spine's
                // heading there turns it a right angle about `y`.
                let want = Vector::new(0.0, along, 0.0);
                assert!(
                    n.dot(want).abs() < 1e-6 && n.magnitude() > 0.5,
                    "{reversed} {corners:?}: {n:?} at {at:?}"
                );
                // `±y × t` with `t` the spine's heading near `at`.
                let t = if at.z <= 5.0 && at.x <= 0.0 {
                    Vector::new(0.0, 0.0, 1.0)
                } else {
                    heading
                };
                let crossed = want.cross(t);
                assert!(
                    n.dot(crossed) > 0.0,
                    "{reversed} {corners:?}: {n:?} at {at:?} against {crossed:?}"
                );
            }
        }
    }
}

/// A spine placed by a location sweeps the sheet the spine where it
/// stands sweeps: the same area, its bound moved with it.
#[test]
fn an_open_sheet_down_a_placed_spline_is_the_sheet_where_it_stands() {
    let leg = s_leg(Point::ORIGIN, Vector::new(0.0, 0.0, 1.0));
    let profile = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0)];
    let mut model = Model::new();
    let spine = chained(&mut model, std::slice::from_ref(&leg));
    let wire = polyline(&mut model, &profile);
    let sheet = sweep(&mut model, &wire, &spine, PipeCorners::Mitre);
    let (want, _) = measured(&model, &sheet);
    let shift = Vector::new(0.0, 0.0, -2.0);
    let datum = model.add_datum(ogeom_math::Transform::translation(-shift));
    let stored = chained(
        &mut model,
        &[ogeom_geom::Transformable::transformed(
            &leg,
            &ogeom_math::Transform::translation(shift),
            T,
        )
        .unwrap()],
    );
    let placed = stored.moved(&ogeom_topo::Location::of(datum));
    let sheet = sweep(&mut model, &wire, &placed, PipeCorners::Mitre);
    let (area, _) = measured(&model, &sheet);
    assert!(
        (area - want).abs() <= want * 1e-12,
        "placed {area} against {want}"
    );
}
