//! Pipe shells under the frame laws beyond the rotation-minimizing and
//! Frenet frames, and the ways of turning a spine's sharp corners.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_edge, make_face, make_polygon, make_wire, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, PlaneSurface};
use ogeom::math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::offset::{PipeCorners, PipeLaw, make_pipe_shell_with};
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

/// A 2 x 2 square centred on the origin in the XZ plane.
fn square(model: &mut Model) -> Shape {
    let pts =
        [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)].map(|(x, z)| Point::new(x, 0.0, z));
    let wire = make_polygon(model, &pts, true, T).unwrap().shape;
    let plane = Plane::new(Frame::new(Point::ORIGIN, -Direction::Y, Direction::X, T).unwrap());
    make_face(model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape
}

/// The quarter arc about (20, 0, 0) of radius 20 in the XY plane, from the
/// origin (heading -Y) to (20, -20, 0).
fn quarter_arc(model: &mut Model) -> Shape {
    let frame = Frame::new(Point::new(20.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let circle = Circle::new(frame, 20.0, T).unwrap();
    let edge = make_edge(model, CircleCurve::new(circle).into(), (PI, 1.5 * PI), T)
        .unwrap()
        .shape;
    make_wire(model, &[edge], T).unwrap().shape
}

fn highest(model: &Model, shape: &Shape, axis: impl Fn(Point) -> f64) -> f64 {
    let mesh =
        ogeom::mesh::triangulate(model, shape, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
    mesh.positions
        .iter()
        .map(|p| axis(*p))
        .fold(f64::NEG_INFINITY, f64::max)
}

#[test]
fn a_leaning_binormal_turns_the_section_as_it_goes() {
    let mut model = Model::new();
    let profile = square(&mut model);
    let spine = quarter_arc(&mut model);
    let law = PipeLaw::Binormal(Direction::new(Vector::new(1.0, 0.0, 1.0), T).unwrap());
    let pipe = make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &law,
        PipeCorners::Mitre,
        1e-3,
        T,
    )
    .unwrap()
    .shape;
    assert!(ogeom::algo::check(&model, &pipe, T).unwrap().is_valid());
    let v = volume(&model, &pipe);
    let want = 4.0 * 10.0 * PI;
    assert!((v - want).abs() < want * 5e-3, "{v} against {want}");
    // Measured at the default deflection too: the walls are splines whose
    // fitted pcurves stray past the chart's border by rounding, and the
    // exact integral takes them all the same.
    let coarse = volume_properties(&model, &pipe, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((coarse - want).abs() < 5e-3, "{coarse} against {want}");
    let top = highest(&model, &pipe, |p| p.z);
    assert!((top - 2.0_f64.sqrt()).abs() < 1e-2, "reaches z {top}");
}

#[test]
fn a_binormal_square_to_a_planar_spine_is_the_rotation_minimizing_pipe() {
    let mut model = Model::new();
    let profile = square(&mut model);
    let spine = quarter_arc(&mut model);
    let pipe = make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::Binormal(Direction::Z),
        PipeCorners::Mitre,
        1e-3,
        T,
    )
    .unwrap()
    .shape;
    let plain = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    let (a, b) = (volume(&model, &pipe), volume(&model, &plain));
    assert!((a - b).abs() < b * 1e-4, "{a} against {b}");
    let top = highest(&model, &pipe, |p| p.z);
    assert!((top - 1.0).abs() < 1e-3, "reaches z {top}");
}

#[test]
fn an_auxiliary_spine_turns_the_section_toward_it() {
    let mut model = Model::new();
    let profile = square(&mut model);
    let spine = make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(0.0, 20.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let guide = make_polygon(
        &mut model,
        &[Point::new(5.0, 0.0, 0.0), Point::new(0.0, 20.0, 5.0)],
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
    assert!(ogeom::algo::check(&model, &pipe, T).unwrap().is_valid());
    let v = volume(&model, &pipe);
    assert!((v - 80.0).abs() < 80.0 * 5e-3, "{v}");
    // Halfway along, the square is turned an eighth of a turn: its corner
    // reaches x = sqrt 2.
    let half = Plane::through(Point::new(0.0, 10.0, 0.0), Direction::Y);
    let cut = ogeom::boolean::section_face(&mut model, &pipe, &half, T)
        .unwrap()
        .shape;
    let reach = highest_on_edges(&model, &cut);
    assert!((reach - 2.0_f64.sqrt()).abs() < 1e-2, "reaches x {reach}");
}

/// A guide ending a few micrometres short of the last station's plane, as
/// single-precision sketch data does, is taken to cross it at its end.
#[test]
fn an_auxiliary_spine_ending_a_hair_short_still_turns_the_section() {
    let mut model = Model::new();
    let profile = square(&mut model);
    let spine = make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(0.0, 20.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let guide = make_polygon(
        &mut model,
        &[
            Point::new(5.0, 0.0, 0.0),
            Point::new(0.0, 20.0 - 3.1e-6, 5.0),
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
    assert!(ogeom::algo::check(&model, &pipe, T).unwrap().is_valid());
    let v = volume(&model, &pipe);
    assert!((v - 80.0).abs() < 80.0 * 5e-3, "{v}");
}

fn highest_on_edges(model: &Model, shape: &Shape) -> f64 {
    use ogeom::geom::Curve3d as _;
    let mut best = f64::NEG_INFINITY;
    for edge in ogeom::topo::explore_unique(model, shape, ogeom::topo::ShapeType::Edge).unwrap() {
        let data = model.node(&edge).unwrap().data().as_edge().unwrap().clone();
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let c = model.geometry().curve(*curve).unwrap();
        for k in 0..=200 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 200.0;
            best = best.max(c.point_at(t, T).unwrap().x);
        }
    }
    best
}

/// The unit disc in the XZ plane about the origin, piped round the L
/// (0,0,0) -> (0,10,0) -> (10,10,0) with each way of turning its corner.
#[test]
fn the_three_ways_round_a_corner() {
    for (corners, want) in [
        (PipeCorners::Mitre, 20.0 * PI),
        (PipeCorners::Extended, 22.0 * PI - 16.0 / 3.0),
        (PipeCorners::Round, 20.0 * PI - 4.0 / 3.0 + PI / 3.0),
    ] {
        let mut model = Model::new();
        let frame = Frame::new(Point::ORIGIN, -Direction::Y, Direction::X, T).unwrap();
        let circle = Circle::new(frame, 1.0, T).unwrap();
        let edge = make_edge(
            &mut model,
            CircleCurve::new(circle).into(),
            (0.0, 2.0 * PI),
            T,
        )
        .unwrap()
        .shape;
        let wire = make_wire(&mut model, &[edge], T).unwrap().shape;
        let disc = make_face(
            &mut model,
            PlaneSurface::new(Plane::new(frame)).into(),
            &[wire],
            T,
        )
        .unwrap()
        .shape;
        let spine = make_polygon(
            &mut model,
            &[
                Point::ORIGIN,
                Point::new(0.0, 10.0, 0.0),
                Point::new(10.0, 10.0, 0.0),
            ],
            false,
            T,
        )
        .unwrap()
        .shape;
        let pipe = make_pipe_shell_with(
            &mut model,
            &disc,
            &spine,
            &PipeLaw::RotationMinimizing,
            corners,
            1e-3,
            T,
        )
        .unwrap_or_else(|e| panic!("{corners:?}: {e}"))
        .shape;
        assert!(
            ogeom::algo::check(&model, &pipe, T).unwrap().is_valid(),
            "{corners:?}"
        );
        let v = volume(&model, &pipe);
        assert!(
            (v - want).abs() < want * 5e-3,
            "{corners:?}: {v} against {want}"
        );
    }
}

fn circle_wire(model: &mut Model, frame: Frame, radius: f64) -> Shape {
    let circle = Circle::new(frame, radius, T).unwrap();
    let edge = make_edge(model, CircleCurve::new(circle).into(), (0.0, 2.0 * PI), T)
        .unwrap()
        .shape;
    make_wire(model, &[edge], T).unwrap().shape
}

/// Two circles down a line, radius 2 at its start and 3 at its end: the
/// frustum between them.
#[test]
fn a_pipe_through_two_sections_changes_its_shape_down_the_path() {
    let mut model = Model::new();
    let small = circle_wire(
        &mut model,
        Frame::new(Point::ORIGIN, Direction::Y, Direction::X, T).unwrap(),
        2.0,
    );
    let large = circle_wire(
        &mut model,
        Frame::new(Point::new(0.0, 20.0, 0.0), Direction::Y, Direction::X, T).unwrap(),
        3.0,
    );
    let spine = make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(0.0, 20.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let pipe =
        ogeom::offset::make_pipe_sections(&mut model, &[small, large], &spine, false, 1e-3, T)
            .unwrap()
            .shape;
    assert!(ogeom::algo::check(&model, &pipe, T).unwrap().is_valid());
    let v = volume(&model, &pipe);
    let want = PI * 20.0 / 3.0 * (4.0 + 6.0 + 9.0);
    assert!((v - want).abs() < want * 1e-2, "{v} against {want}");
}

/// The end circle drawn with its seam on the other side, turned a quarter,
/// or running the other way, as a face or as a wire: the sections' starts
/// are matched before blending, so every one is the same frustum.
#[test]
fn sections_drawn_from_different_starts_blend_without_a_twist() {
    let want = PI * 20.0 / 3.0 * (4.0 + 6.0 + 9.0);
    let ends = [
        (Direction::Y, -Direction::X, false),
        (Direction::Y, -Direction::X, true),
        (Direction::Y, Direction::Z, false),
        (-Direction::Y, Direction::X, false),
        (-Direction::Y, -Direction::Z, true),
    ];
    for (normal, x_axis, as_face) in ends {
        let mut model = Model::new();
        let small = circle_wire(
            &mut model,
            Frame::new(Point::ORIGIN, Direction::Y, Direction::X, T).unwrap(),
            2.0,
        );
        let frame = Frame::new(Point::new(0.0, 20.0, 0.0), normal, x_axis, T).unwrap();
        let mut large = circle_wire(&mut model, frame, 3.0);
        let mut small = small;
        if as_face {
            let face = |model: &mut Model, wire: Shape, frame: Frame| {
                make_face(
                    model,
                    PlaneSurface::new(Plane::new(frame)).into(),
                    &[wire],
                    T,
                )
                .unwrap()
                .shape
            };
            let start = Frame::new(Point::ORIGIN, Direction::Y, Direction::X, T).unwrap();
            small = face(&mut model, small, start);
            large = face(&mut model, large, frame);
        }
        let spine = make_polygon(
            &mut model,
            &[Point::ORIGIN, Point::new(0.0, 20.0, 0.0)],
            false,
            T,
        )
        .unwrap()
        .shape;
        let pipe =
            ogeom::offset::make_pipe_sections(&mut model, &[small, large], &spine, false, 1e-3, T)
                .unwrap()
                .shape;
        assert!(ogeom::algo::check(&model, &pipe, T).unwrap().is_valid());
        let v = volume(&model, &pipe);
        assert!(
            (v - want).abs() < want * 1e-2,
            "{v} against {want} for {normal:?}, {x_axis:?}"
        );
    }
}

/// Two equal circles at the ends of a quarter arc: the same tube a single
/// section sweeps, its area times the arc.
#[test]
fn two_equal_sections_round_an_arc_are_the_plain_pipe() {
    let mut model = Model::new();
    let spine = quarter_arc(&mut model);
    let first = circle_wire(
        &mut model,
        Frame::new(Point::ORIGIN, -Direction::Y, Direction::X, T).unwrap(),
        2.0,
    );
    let last = circle_wire(
        &mut model,
        Frame::new(Point::new(20.0, -20.0, 0.0), Direction::X, Direction::Y, T).unwrap(),
        2.0,
    );
    let pipe =
        ogeom::offset::make_pipe_sections(&mut model, &[first, last], &spine, false, 1e-3, T)
            .unwrap()
            .shape;
    assert!(ogeom::algo::check(&model, &pipe, T).unwrap().is_valid());
    let v = volume(&model, &pipe);
    let want = PI * 4.0 * 10.0 * PI;
    assert!((v - want).abs() < want * 5e-3, "{v} against {want}");
}

/// A square carried along a quarter circle without turning: every level
/// cut is the band a 2-long segment sweeps along the arc, so the solid is
/// 2 x 20 x 2, and the section at the end still stands square to Y, where
/// the arc runs along X. Along a straight spine it is the plain pipe.
#[test]
fn a_fixed_section_is_carried_without_turning() {
    let mut model = Model::new();
    let profile = square(&mut model);
    // About (20, 0, 0), from the origin heading +Y round to (20, 20, 0).
    let frame = Frame::new(Point::new(20.0, 0.0, 0.0), -Direction::Z, -Direction::X, T).unwrap();
    let circle = Circle::new(frame, 20.0, T).unwrap();
    let edge = make_edge(
        &mut model,
        CircleCurve::new(circle).into(),
        (0.0, PI / 2.0),
        T,
    )
    .unwrap()
    .shape;
    let spine = make_wire(&mut model, &[edge], T).unwrap().shape;
    let pipe = make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::Fixed,
        PipeCorners::Mitre,
        1e-3,
        T,
    )
    .unwrap()
    .shape;
    let diagnosis = ogeom::algo::check(&model, &pipe, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let v = volume(&model, &pipe);
    assert!((v - 80.0).abs() < 80.0 * 1e-3, "{v}");
    // The far end: the section at y = 20, from x = 19 to 21, not turned.
    let mesh =
        ogeom::mesh::triangulate(&model, &pipe, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
    let far: Vec<Point> = mesh
        .positions
        .iter()
        .copied()
        .filter(|p| (p.y - 20.0).abs() < 1e-9)
        .collect();
    let (lo, hi) = far
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
            (lo.min(p.x), hi.max(p.x))
        });
    assert!(
        (lo - 19.0).abs() < 1e-6 && (hi - 21.0).abs() < 1e-6,
        "{lo} .. {hi}"
    );
    let top = mesh
        .positions
        .iter()
        .map(|p| p.y)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(top < 20.0 + 1e-9, "reaches y {top}");

    // Straight, it is the plain pipe: the same prism.
    let line = make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(0.0, 10.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let (fixed, plain) = (
        make_pipe_shell_with(
            &mut model,
            &profile,
            &line,
            &PipeLaw::Fixed,
            PipeCorners::Mitre,
            1e-3,
            T,
        )
        .unwrap()
        .shape,
        ogeom::offset::make_pipe_shell(&mut model, &profile, &line, false, 1e-3, T)
            .unwrap()
            .shape,
    );
    assert!(ogeom::algo::check(&model, &fixed, T).unwrap().is_valid());
    let (a, b) = (volume(&model, &fixed), volume(&model, &plain));
    assert!(
        (a - b).abs() < 1e-9 && (a - 40.0).abs() < 1e-9,
        "{a} against {b}"
    );
}

/// A one by a quarter rectangle up to 2e7 from the origin on a tilted
/// plane, carried without turning along a straight spine of ten square to
/// it: the section's normal is its own ring's vector area, so the solid is
/// the prism there as at the origin.
#[test]
fn a_fixed_section_far_from_the_origin_sweeps_its_prism() {
    let n = Direction::new(Vector::new(1.0, 1.0, 1.0), T).unwrap();
    for o in [0.0, 1e6, 1e7, 2e7] {
        let mut model = Model::new();
        let frame = Frame::new(Point::new(o, -o, o), n, Direction::X, T).unwrap();
        let at = |x: f64, y: f64| frame.origin() + frame.x().vector() * x + frame.y().vector() * y;
        let ring = [at(0.0, 0.0), at(1.0, 0.0), at(1.0, 0.25), at(0.0, 0.25)];
        let wire = make_polygon(&mut model, &ring, true, T).unwrap().shape;
        let profile = make_face(
            &mut model,
            PlaneSurface::new(Plane::new(frame)).into(),
            &[wire],
            T,
        )
        .unwrap()
        .shape;
        let spine = make_polygon(
            &mut model,
            &[at(0.5, 0.125), at(0.5, 0.125) + frame.z().vector() * 10.0],
            false,
            T,
        )
        .unwrap()
        .shape;
        let pipe = make_pipe_shell_with(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::Fixed,
            PipeCorners::Mitre,
            1e-3,
            T,
        )
        .unwrap_or_else(|e| panic!("o {o:e}: {e}"))
        .shape;
        let diagnosis = ogeom::algo::check(&model, &pipe, T).unwrap();
        assert!(diagnosis.is_valid(), "o {o:e}: {diagnosis}");
        let v = volume_properties(&model, &pipe, Deflection::default(), T)
            .unwrap_or_else(|e| panic!("o {o:e}: {e}"))
            .mass;
        assert!((v - 2.5).abs() < 2.5 * 1e-6, "o {o:e}: volume {v}");
    }
}
