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

/// The 2 x 2 square turning a corner of next to nothing extended, as far
/// down as an angle's tolerance: each leg's prism mitres to `4 l`, and the
/// squared outside of the corner adds `2 tan(θ/2)`. Every point is held to
/// the confusion distance, so the volume to the solid's area times it.
#[test]
fn a_square_turning_next_to_nothing_extends_to_its_legs_volume() {
    let l = 10.0;
    for theta in [1e-11_f64, 1e-9, 1e-7] {
        let mut model = Model::new();
        let face = square(&mut model);
        let spine = make_polygon(
            &mut model,
            &[
                Point::ORIGIN,
                Point::new(0.0, l, 0.0),
                Point::new(l * theta.sin(), l + l * theta.cos(), 0.0),
            ],
            false,
            T,
        )
        .unwrap()
        .shape;
        let pipe = make_pipe_shell_with(
            &mut model,
            &face,
            &spine,
            &PipeLaw::RotationMinimizing,
            PipeCorners::Extended,
            1e-3,
            T,
        )
        .unwrap_or_else(|e| panic!("{theta}: {e}"))
        .shape;
        assert!(
            ogeom::algo::check(&model, &pipe, T).unwrap().is_valid(),
            "{theta}"
        );
        let want = 8.0 * l + 2.0 * (theta / 2.0).tan();
        let area = ogeom::algo::surface_properties(
            &model,
            &pipe,
            Deflection::with_chord(1e-3).unwrap(),
            T,
        )
        .unwrap()
        .mass;
        let v = volume(&model, &pipe);
        assert!(
            (v - want).abs() <= area * T.confusion(),
            "{theta}: {v} against {want}"
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

/// A square of side `2 half` with a circular hole of `radius`, both about
/// `(0, y, 0)` in the plane square to Y there.
fn holed_square(model: &mut Model, y: f64, half: f64, radius: f64) -> Shape {
    let pts = [(-half, -half), (half, -half), (half, half), (-half, half)]
        .map(|(x, z)| Point::new(x, y, z));
    let outer = make_polygon(model, &pts, true, T).unwrap().shape;
    let frame = Frame::new(Point::new(0.0, y, 0.0), Direction::Y, Direction::X, T).unwrap();
    let hole = circle_wire(model, frame, radius);
    let plane = Plane::new(frame);
    make_face(model, PlaneSurface::new(plane).into(), &[outer, hole], T)
        .unwrap()
        .shape
}

/// Two holed squares down a line, the second half the size: the outer
/// loops blend into the pyramid's frustum and the holes into the cone's,
/// cut from it. Each solid's volume is `h (A₀ + √(A₀ A₁) + A₁) / 3`; the
/// measured volume strays from the difference by at most the skin's
/// tolerance and the mesh's chord over the walls' area, bounded by each
/// loop's larger perimeter times the slant run.
#[test]
fn a_pipe_through_holed_sections_keeps_its_hole_down_the_path() {
    let (h, skin, chord) = (20.0, 1e-3, 1e-3);
    let mut model = Model::new();
    let first = holed_square(&mut model, 0.0, 2.0, 1.0);
    let last = holed_square(&mut model, h, 1.0, 0.5);
    let spine = make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(0.0, h, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let pipe =
        ogeom::offset::make_pipe_sections(&mut model, &[first, last], &spine, false, skin, T)
            .unwrap()
            .shape;
    let found = ogeom::algo::check(&model, &pipe, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let frustum = |a0: f64, a1: f64| h * (a0 + (a0 * a1).sqrt() + a1) / 3.0;
    let want = frustum(16.0, 4.0) - frustum(PI, PI * 0.25);
    let area = 16.0 * h.hypot(1.0) + 2.0 * PI * h.hypot(0.5);
    let bound = (skin + chord) * area;
    let v = volume(&model, &pipe);
    assert!(
        (v - want).abs() < bound,
        "{v} against {want}, bound {bound}"
    );
}

/// A holed square round the quarter arc under the leaning binormal: the
/// hole sweeps its own solid under the same frames and is cut from the
/// outer one. Both centroids ride the spine, so Pappus gives the holed
/// area times the arc's length; the measured volume strays from it by at
/// most the skin's tolerance and the mesh's chord over each wall's area,
/// its loop's perimeter times its longest path round the arc.
#[test]
fn a_holed_square_under_a_leaning_binormal_keeps_its_hole() {
    let (half, radius, big, skin, chord) = (1.0, 0.5, 20.0, 1e-3, 1e-3);
    let mut model = Model::new();
    let profile = holed_square(&mut model, 0.0, half, radius);
    let spine = quarter_arc(&mut model);
    let law = PipeLaw::Binormal(Direction::new(Vector::new(1.0, 0.0, 1.0), T).unwrap());
    let pipe = make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &law,
        PipeCorners::Mitre,
        skin,
        T,
    )
    .unwrap()
    .shape;
    let found = ogeom::algo::check(&model, &pipe, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let length = big * PI / 2.0;
    let want = (4.0 * half * half - PI * radius * radius) * length;
    let path = |reach: f64| length * (1.0 + reach / big);
    let area = 8.0 * half * path(half * 2.0_f64.sqrt()) + 2.0 * PI * radius * path(radius);
    let bound = (skin + chord) * area;
    let v = volume(&model, &pipe);
    assert!(
        (v - want).abs() < bound,
        "{v} against {want}, bound {bound}"
    );
}

#[test]
fn a_pipe_through_sections_with_different_hole_counts_is_refused_by_name() {
    let mut model = Model::new();
    let first = holed_square(&mut model, 0.0, 2.0, 1.0);
    let last = square(&mut model);
    let spine = make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(0.0, -20.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let err = ogeom::offset::make_pipe_sections(&mut model, &[first, last], &spine, false, 1e-3, T)
        .unwrap_err();
    assert!(err.to_string().contains("as many holes"), "{err}");
}

/// The spine in the `xy` plane: the quarter arc of radius `big` about
/// `(-big, 0, 0)` from the origin, heading `+y`, to `(-big, big, 0)`,
/// heading `-x`; then a sharp quarter turn up the line to
/// `(-big, big + run, 0)`.
fn arc_then_turn(model: &mut Model, big: f64, run: f64) -> Shape {
    arc_then_turn_at(model, big, run, 0.0)
}

/// [`arc_then_turn`] raised to the height `z`.
fn arc_then_turn_at(model: &mut Model, big: f64, run: f64, z: f64) -> Shape {
    use ogeom::geom::Curve3d as _;
    let frame = Frame::new(Point::new(-big, 0.0, z), Direction::Z, Direction::X, T).unwrap();
    let circle = Circle::new(frame, big, T).unwrap();
    let arc = make_edge(model, CircleCurve::new(circle).into(), (0.0, 0.5 * PI), T)
        .unwrap()
        .shape;
    let (_, joint) = ogeom::algo::edge_vertices(model, &arc).unwrap().unwrap();
    let (from, to) = (Point::new(-big, big, z), Point::new(-big, big + run, z));
    let end = ogeom::algo::make_vertex(model, to).shape;
    let segment = ogeom::geom::LineCurve::segment(from, to, T).unwrap();
    let range = segment.domain();
    let line = ogeom::algo::make_edge_between(model, segment.into(), range, &joint, &end, T)
        .unwrap()
        .shape;
    make_wire(model, &[arc, line], T).unwrap().shape
}

fn face_count(model: &Model, shape: &Shape) -> usize {
    ogeom::topo::explore_unique(model, shape, ogeom::topo::ShapeType::Face)
        .unwrap()
        .len()
}

/// A holed square down an arc into a sharp quarter turn, in the `xy`
/// plane, under the binormal `z`: on a planar spine that law's frames are
/// the rotation-minimizing ones, so each corner it builds is the
/// rotation-minimizing sweep's. The arc's run is skinned under the law
/// where the rotation-minimizing frame sweeps it exactly; every other
/// piece (the line's prism, the corner's extensions, the holes cut) is
/// built alike from the same sections. So the two solids differ by at
/// most the skin's tolerance over the arc run's walls, and each strays
/// from its own integral by at most the mesh's chord over its walls: the
/// bound is their sum over a bound on the walls' area, the perimeters
/// times the longest path a profile point runs, round the arc, up the
/// line and past the corner by the section's reach either way.
///
/// The extended corner matches in volume, face count and validity. The
/// round corner matches in volume and validity: the next run starts on
/// the arc run's end section turned within the skin's tolerance, which
/// is what the skinned run is held to.
#[test]
fn a_holed_square_under_a_level_binormal_turns_a_corner_as_the_rotation_minimizing_sweep() {
    let (big, run, half, radius, skin, chord) = (20.0, 10.0, 1.0, 0.5, 1e-3, 1e-3);
    let sweep = |law: &PipeLaw<'_>, corners: PipeCorners| {
        let mut model = Model::new();
        let profile = holed_square(&mut model, 0.0, half, radius);
        let spine = arc_then_turn(&mut model, big, run);
        let built = make_pipe_shell_with(&mut model, &profile, &spine, law, corners, skin, T)
            .map(|b| b.shape);
        (model, built)
    };
    let level = PipeLaw::Binormal(Direction::Z);

    let (rmf_model, rmf) = sweep(&PipeLaw::RotationMinimizing, PipeCorners::Extended);
    let rmf = rmf.unwrap();
    let (law_model, swept) = sweep(&level, PipeCorners::Extended);
    let swept = swept.unwrap();
    let found = ogeom::algo::check(&law_model, &swept, T).unwrap();
    assert!(found.is_valid(), "{found}");
    assert!(ogeom::algo::check(&rmf_model, &rmf, T).unwrap().is_valid());
    assert_eq!(face_count(&law_model, &swept), face_count(&rmf_model, &rmf));
    let reach = half * 2.0_f64.sqrt();
    let path = big * PI / 2.0 * (1.0 + reach / big) + run + 2.0 * reach;
    let perimeter = 8.0 * half + 2.0 * PI * radius;
    let bound = (skin + 2.0 * chord) * perimeter * path;
    let (v, want) = (volume(&law_model, &swept), volume(&rmf_model, &rmf));
    assert!(
        (v - want).abs() < bound,
        "{v} against {want}, bound {bound}"
    );

    let (rmf_model, rmf) = sweep(&PipeLaw::RotationMinimizing, PipeCorners::Round);
    let rmf = rmf.unwrap();
    let (law_model, swept) = sweep(&level, PipeCorners::Round);
    let swept = swept.unwrap();
    let found = ogeom::algo::check(&law_model, &swept, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let (v, want) = (volume(&law_model, &swept), volume(&rmf_model, &rmf));
    assert!(
        (v - want).abs() < bound,
        "round: {v} against {want}, bound {bound}"
    );
}

/// The spine of [`arc_then_turn`] walked the other way: up the line from
/// `(-big, big + run, 0)` heading `-y` to `(-big, big, 0)`, then a sharp
/// quarter turn onto the quarter arc, heading `+x`, back to the origin.
fn turn_then_arc(model: &mut Model, big: f64, run: f64) -> Shape {
    use ogeom::geom::Curve3d as _;
    let frame = Frame::new(Point::new(-big, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let circle = Circle::new(frame, big, T).unwrap();
    let arc = make_edge(model, CircleCurve::new(circle).into(), (0.0, 0.5 * PI), T)
        .unwrap()
        .shape;
    let (_, joint) = ogeom::algo::edge_vertices(model, &arc).unwrap().unwrap();
    let (from, to) = (Point::new(-big, big + run, 0.0), Point::new(-big, big, 0.0));
    let start = ogeom::algo::make_vertex(model, from).shape;
    let segment = ogeom::geom::LineCurve::segment(from, to, T).unwrap();
    let range = segment.domain();
    let line = ogeom::algo::make_edge_between(model, segment.into(), range, &start, &joint, T)
        .unwrap()
        .shape;
    make_wire(model, &[line, arc.reversed()], T).unwrap().shape
}

/// A profile square to `-y` at `(-big, big + run, 0)`, where
/// [`turn_then_arc`] starts: the holed square, or the disc of radius
/// `half` where `radius` is naught.
fn profile_down_y(model: &mut Model, at: Point, half: f64, radius: f64) -> Shape {
    let frame = Frame::new(at, -Direction::Y, Direction::X, T).unwrap();
    let plane = PlaneSurface::new(Plane::new(frame));
    if radius == 0.0 {
        let disc = circle_wire(model, frame, half);
        return make_face(model, plane.into(), &[disc], T).unwrap().shape;
    }
    let (u, v) = (frame.x().vector(), frame.y().vector());
    let pts = [(-half, -half), (half, -half), (half, half), (-half, half)]
        .map(|(a, b)| at + u * a + v * b);
    let outer = make_polygon(model, &pts, true, T).unwrap().shape;
    let hole = circle_wire(model, frame, radius);
    make_face(model, plane.into(), &[outer, hole], T)
        .unwrap()
        .shape
}

/// The disc of radius `half` about `(0, y, 0)` in the plane square to Y.
fn disc(model: &mut Model, y: f64, half: f64) -> Shape {
    let frame = Frame::new(Point::new(0.0, y, 0.0), Direction::Y, Direction::X, T).unwrap();
    let wire = circle_wire(model, frame, half);
    make_face(
        model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape
}

/// The mitred solid's volume where an arc of radius `big` meets a line
/// `run` long at a right angle, the corner turning against the arc, for a
/// section of area `area` centred on the spine whose height across the
/// turn's plane at a distance `u` either side of the spine is
/// `width(u)`, out to `reach`.
///
/// Away from the corner each run holds its area times its length
/// (Pappus round the arc). The line's run is run back past the corner
/// and cut on the plane halving it, which takes inside the corner what it
/// adds outside, the section being symmetric. The arc's run is run on by
/// its end section's prism, which adds a wedge `u` deep at `u` outside
/// the corner; inside, the plane cuts the arc's own run: at `u` there, on
/// the circle of radius `ρ = big + u`, the plane through the corner at
/// half a right angle to the end section meets that circle `ψ(ρ)` back
/// round the arc, where with the arc's centre at the origin and the
/// corner at `(big, 0)` its point `(big + t, t)` has
/// `2t² + 2 big t + big² - ρ² = 0`, so `ρ ψ` of arc is cut. Each column
/// `u` adds `(u - ρ ψ(ρ)) width(u)`, integrated by the midpoint rule over
/// a hundred thousand columns.
fn mitred_against_an_arc(
    big: f64,
    run: f64,
    area: f64,
    reach: f64,
    width: impl Fn(f64) -> f64,
) -> f64 {
    const COLUMNS: u32 = 100_000;
    let psi = |rho: f64| {
        let t = (-big + 2.0f64.mul_add(rho * rho, -big * big).sqrt()) / 2.0;
        t.atan2(big + t)
    };
    let du = reach / f64::from(COLUMNS);
    let corner: f64 = (0..COLUMNS)
        .map(|k| {
            let u = (f64::from(k) + 0.5) * du;
            let rho = big + u;
            (u - rho * psi(rho)) * width(u) * du
        })
        .sum();
    area * (big * PI / 2.0 + run) + corner
}

/// A section with a hole, or a disc, mitred where a skinned run meets a
/// straight one, under the binormal and auxiliary laws: the arc's run is
/// skinned and ends on the skin's fit of each loop, which wobbles across
/// the exact loop the law places on the line's run. The line's prism takes
/// the skin's own fit of each loop the two share (the hole, the disc;
/// not the square, which the leaning binormal turns there), so both cut
/// the mitre plane alike and the pieces fuse. Every section here is
/// level at the corner and symmetric across the spine, so each solid is
/// [`mitred_against_an_arc`]'s.
///
/// Arc first, the arc's run is the arc swept alone, measured for its own
/// skin's error; what is left strays by at most the skin's tolerance over
/// the walls the corner's wedges run along and over the line's walls,
/// whose loops are fits. Line first, the arc's run is skinned from the
/// other end, so its whole wall counts as well.
#[test]
fn a_smooth_loop_mitres_against_a_skinned_run_under_a_law() {
    let (big, run, half, radius, skin) = (20.0, 10.0, 1.0, 0.5, 1e-3);
    let reach = half * 2.0_f64.sqrt();
    let holed = (
        4.0 * half * half - PI * radius * radius,
        8.0 * half + 2.0 * PI * radius,
    );
    let round = (PI * half * half, 2.0 * PI * half);
    let holed_width = |u: f64| {
        if u > half {
            0.0
        } else if u < radius {
            2.0 * half - 2.0 * (radius * radius - u * u).sqrt()
        } else {
            2.0 * half
        }
    };
    let disc_width = |u: f64| 2.0 * (half * half - u * u).max(0.0).sqrt();
    let lean = Direction::new(Vector::new(1.0, 0.0, 1.0), T).unwrap();
    for (case, has_hole, arc_first) in [
        ("holed, level", true, true),
        ("holed, leaning", true, true),
        ("disc, leaning", false, true),
        ("holed, auxiliary", true, true),
        ("holed, level, line first", true, false),
    ] {
        let mut model = Model::new();
        let (profile, spine) = if arc_first {
            let profile = if has_hole {
                holed_square(&mut model, 0.0, half, radius)
            } else {
                disc(&mut model, 0.0, half)
            };
            (profile, arc_then_turn(&mut model, big, run))
        } else {
            let at = Point::new(-big, big + run, 0.0);
            let profile = profile_down_y(&mut model, at, half, radius);
            (profile, turn_then_arc(&mut model, big, run))
        };
        let guide = arc_then_turn_at(&mut model, big, run, 5.0);
        let law = if case.contains("leaning") {
            PipeLaw::Binormal(lean)
        } else if case.contains("auxiliary") {
            PipeLaw::Auxiliary { guide: &guide }
        } else {
            PipeLaw::Binormal(Direction::Z)
        };
        let pipe = make_pipe_shell_with(
            &mut model,
            &profile,
            &spine,
            &law,
            PipeCorners::Mitre,
            skin,
            T,
        )
        .unwrap_or_else(|e| panic!("{case}: {e}"))
        .shape;
        let found = ogeom::algo::check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{case}: {found}");
        let v = volume(&model, &pipe);
        let (area, perimeter) = if has_hole { holed } else { round };
        let want = if has_hole {
            mitred_against_an_arc(big, run, area, reach, holed_width)
        } else {
            mitred_against_an_arc(big, run, area, half, disc_width)
        };
        let corner_walls = perimeter * (2.0 * reach + run);
        let (own, bound) = if arc_first {
            let arc = model.ordered_children_of(&spine).unwrap()[0].clone();
            let arc = make_wire(&mut model, &[arc], T).unwrap().shape;
            let alone = make_pipe_shell_with(
                &mut model,
                &profile,
                &arc,
                &law,
                PipeCorners::Mitre,
                skin,
                T,
            )
            .unwrap()
            .shape;
            let own = volume(&model, &alone) - area * big * PI / 2.0;
            (own, skin * corner_walls)
        } else {
            let arc_walls = perimeter * (big + reach) * PI / 2.0;
            (0.0, skin * (corner_walls + arc_walls))
        };
        assert!(
            (v - own - want).abs() < bound,
            "{case}: {v} less the arc's own {own:e} against {want}, bound {bound}"
        );
    }
}

/// A closed square wire, `2a` across, up a leg `l` long and a quarter
/// turn along another, turned run by run: the walls of the square face's
/// solid, the end caps taken off. Extended or mitred, each side wall in
/// the turn's plane is the L the legs make, `4 a l`, the outer walls run
/// `l + a` and the inner `l - a`, each `2a` wide: `16 a l` in all.
/// Rounded, the legs end square at the corner: each side wall loses the
/// outside corner's `a²` and gains a quarter disc `π a² / 4`, the outer
/// walls run `l` each and the quarter drum `π a²` joins them:
/// `16 a l - 6 a² + 3π a² / 2`. Under the binormal square to the legs'
/// plane the sections are those the rotation-minimizing frame carries.
#[test]
fn a_closed_wire_turns_corners_run_by_run_as_the_walls_of_its_face() {
    let (a, l) = (1.0, 10.0);
    let squared = 16.0 * a * l;
    let rounded = squared - 6.0 * a * a + 1.5 * PI * a * a;
    for (law, corners, want) in [
        (PipeLaw::RotationMinimizing, PipeCorners::Extended, squared),
        (PipeLaw::RotationMinimizing, PipeCorners::Round, rounded),
        (
            PipeLaw::Binormal(Direction::Y),
            PipeCorners::Extended,
            squared,
        ),
        (PipeLaw::Binormal(Direction::Y), PipeCorners::Round, rounded),
        (PipeLaw::Binormal(Direction::Y), PipeCorners::Mitre, squared),
    ] {
        let mut model = Model::new();
        let pts = [(-a, -a), (a, -a), (a, a), (-a, a)].map(|(x, y)| Point::new(x, y, 0.0));
        let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
        let spine = make_polygon(
            &mut model,
            &[
                Point::ORIGIN,
                Point::new(0.0, 0.0, l),
                Point::new(l, 0.0, l),
            ],
            false,
            T,
        )
        .unwrap()
        .shape;
        let walls = make_pipe_shell_with(&mut model, &wire, &spine, &law, corners, 1e-3, T)
            .unwrap_or_else(|e| panic!("{law:?} {corners:?}: {e}"))
            .shape;
        assert_eq!(
            model.kind_of(&walls).unwrap(),
            ogeom::topo::ShapeType::Shell,
            "{law:?} {corners:?}"
        );
        let area = ogeom::algo::surface_properties(
            &model,
            &walls,
            Deflection::with_chord(1e-4).unwrap(),
            T,
        )
        .unwrap();
        assert!(area.deflection == 0.0, "{law:?} {corners:?}: not exact");
        assert!(
            (area.mass - want).abs() < want * 1e-12,
            "{law:?} {corners:?}: {} against {want}",
            area.mass
        );
    }
}

/// The quarter arc of radius `big` about `(x + big, y, 0)` from `(x, y, 0)`,
/// heading `+y`, clockwise about `z` to `(x + big, y + big, 0)`; it starts
/// on `from` where given.
fn clockwise_arc_from(model: &mut Model, x: f64, y: f64, big: f64, from: Option<&Shape>) -> Shape {
    let frame = Frame::new(Point::new(x + big, y, 0.0), Direction::Z, Direction::X, T).unwrap();
    let circle = CircleCurve::new(Circle::new(frame, big, T).unwrap());
    let to = ogeom::algo::make_vertex(model, Point::new(x + big, y + big, 0.0)).shape;
    let from = from
        .cloned()
        .unwrap_or_else(|| ogeom::algo::make_vertex(model, Point::new(x, y, 0.0)).shape);
    ogeom::algo::make_edge_between(model, circle.into(), (0.5 * PI, PI), &to, &from, T)
        .unwrap()
        .shape
        .reversed()
}

/// The spine in the `xy` plane: [`arc_then_turn`]'s quarter arc of radius
/// `big`, from the origin heading `+y` round to `(-big, big, 0)` heading
/// `-x`, then a sharp quarter turn onto [`clockwise_arc_from`] there, of
/// radius `small`, heading `+y`. The first arc turns against the corner,
/// the second with it.
fn arc_then_arc(model: &mut Model, big: f64, small: f64) -> Shape {
    let frame = Frame::new(Point::new(-big, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let circle = Circle::new(frame, big, T).unwrap();
    let first = make_edge(model, CircleCurve::new(circle).into(), (0.0, 0.5 * PI), T)
        .unwrap()
        .shape;
    let (_, joint) = ogeom::algo::edge_vertices(model, &first).unwrap().unwrap();
    let second = clockwise_arc_from(model, -big, big, small, Some(&joint));
    make_wire(model, &[first, second], T).unwrap().shape
}

/// How much of an arc's run of radius `big` the plane halving a corner
/// cuts off, along the circle of radius `rho` a column inside the corner
/// runs round, where the corner turns through twice `half`: with the arc's
/// centre at the origin, the corner at `(big, 0)` and the run ending
/// heading `y`, the plane meets that circle at `(big ± w cos half, -w sin
/// half)`, `+` where the corner turns against the arc (inside the corner
/// is away from its centre), `-` where it turns with it, so
/// `w² ± 2 big w cos half + big² - ρ² = 0`, and `ρ ψ` is cut where `ψ` is
/// the angle back round to that point.
fn cut_round_an_arc(big: f64, rho: f64, half: f64, against: bool) -> f64 {
    let sign = if against { 1.0 } else { -1.0 };
    let (sin, cos) = half.sin_cos();
    let w = sign * ((rho * rho - big * big * sin * sin).sqrt() - big * cos);
    rho * (w * sin).atan2((sign * w).mul_add(cos, big))
}

/// The mitred solid's volume where an arc of radius `big`, turning
/// against the corner, meets an arc of radius `small` turning with it, at
/// a right angle, for sections of area `area` centred on the spine, each
/// symmetric across it, whose heights across the turn's plane at a
/// distance `u` either side of the spine are `widths.0(u)` at the end of
/// the first arc and `widths.1(u)` at the start of the second, out to
/// `reach`.
///
/// Away from the corner each run holds its area times its length
/// (Pappus). Each run is run on past the corner by its end section's
/// prism, which adds a wedge `u` deep at `u` outside the corner, and cut
/// by the plane halving the corner, which inside it takes off the
/// [`cut_round_an_arc`] of that run's own walls: `ρ = big + u` for the
/// first, `ρ = small - u` for the second. Each column `u` adds
/// `(u - cut₁) width₁(u) + (u - cut₂) width₂(u)`, integrated by the
/// midpoint rule over a hundred thousand columns.
fn mitred_between_arcs(
    big: f64,
    small: f64,
    area: f64,
    reach: f64,
    widths: (&dyn Fn(f64) -> f64, &dyn Fn(f64) -> f64),
) -> f64 {
    const COLUMNS: u32 = 100_000;
    let du = reach / f64::from(COLUMNS);
    let corner: f64 = (0..COLUMNS)
        .map(|k| {
            let u = (f64::from(k) + 0.5) * du;
            let first = (u - cut_round_an_arc(big, big + u, PI / 4.0, true)) * widths.0(u);
            let second = (u - cut_round_an_arc(small, small - u, PI / 4.0, false)) * widths.1(u);
            (first + second) * du
        })
        .sum();
    area * (big + small) * PI / 2.0 + corner
}

/// [`arc_then_turn`] with an auxiliary spine that turns the section a
/// quarter turn up the line: [`arc_then_turn_at`] raised by `rise` to the
/// arc's end, then straight to `rise` behind the line's end, so the law's
/// normal turns from `z` to `-x` along the line. Behind the line, the
/// guide crosses no plane square to the arc.
fn guide_turning_the_line(model: &mut Model, big: f64, run: f64, rise: f64) -> Shape {
    let frame = Frame::new(Point::new(-big, 0.0, rise), Direction::Z, Direction::X, T).unwrap();
    let circle = Circle::new(frame, big, T).unwrap();
    let arc = make_edge(model, CircleCurve::new(circle).into(), (0.0, 0.5 * PI), T)
        .unwrap()
        .shape;
    let (_, joint) = ogeom::algo::edge_vertices(model, &arc).unwrap().unwrap();
    let to = Point::new(-big - rise, big + run, 0.0);
    let end = ogeom::algo::make_vertex(model, to).shape;
    let segment = ogeom::geom::LineCurve::segment(Point::new(-big, big, rise), to, T).unwrap();
    let range = ogeom::geom::Curve3d::domain(&segment);
    let line = ogeom::algo::make_edge_between(model, segment.into(), range, &joint, &end, T)
        .unwrap()
        .shape;
    make_wire(model, &[arc, line], T).unwrap().shape
}

/// Two skinned runs mitred where both carry one loop of the section, under
/// the binormal and auxiliary laws: an arc into an arc, a disc and a holed
/// square under the level and the leaning binormal, and an arc into a line
/// the auxiliary spine turns, so the line is skinned too. Each run's skin
/// ends on the exact section the law places there, the circle in its
/// rational form, and the second run starts on the first's end circle
/// mirrored across the mitre plane: the two runs' prisms cut that plane
/// in one ellipse. The square's sides are straight in either skin.
///
/// Every section here is centred on the spine and symmetric across it. At
/// the corner the first run's is level; under the leaning binormal the
/// second arc starts with the square turned an eighth of a turn, and the
/// disc any way round. So each arc into an arc is
/// [`mitred_between_arcs`]'s, and the arc into the line
/// [`mitred_against_an_arc`]'s.
///
/// The volume is integrated on the faces exactly, so it strays from the
/// closed form by at most the skin's tolerance over the skinned walls: the
/// perimeter times the longest path a profile point runs, round each arc
/// at its radius and the section's reach, up the line, and past the
/// corner by the reach either way on each run.
#[test]
fn two_skinned_runs_mitre_on_their_exact_common_section() {
    let (big, small, run, half, radius, skin) = (20.0, 15.0, 10.0, 1.0, 0.5, 1e-3);
    let reach = half * 2.0_f64.sqrt();
    let holed = (
        4.0 * half * half - PI * radius * radius,
        8.0 * half + 2.0 * PI * radius,
    );
    let round = (PI * half * half, 2.0 * PI * half);
    let hole = |u: f64| 2.0 * (radius * radius - u * u).max(0.0).sqrt();
    let holed_width = |u: f64| if u > half { 0.0 } else { 2.0 * half - hole(u) };
    let diamond_width = |u: f64| 2.0 * (reach - u).max(0.0) - hole(u);
    let disc_width = |u: f64| 2.0 * (half * half - u * u).max(0.0).sqrt();
    let lean = Direction::new(Vector::new(1.0, 0.0, 1.0), T).unwrap();
    for (case, has_hole) in [
        ("disc, level", false),
        ("holed, level", true),
        ("disc, leaning", false),
        ("holed, leaning", true),
        ("disc, auxiliary turning the line", false),
    ] {
        let line = case.contains("line");
        let mut model = Model::new();
        let guide = guide_turning_the_line(&mut model, big, run, 5.0);
        let law = if case.contains("leaning") {
            PipeLaw::Binormal(lean)
        } else if line {
            PipeLaw::Auxiliary { guide: &guide }
        } else {
            PipeLaw::Binormal(Direction::Z)
        };
        let profile = if has_hole {
            holed_square(&mut model, 0.0, half, radius)
        } else {
            disc(&mut model, 0.0, half)
        };
        let spine = if line {
            arc_then_turn(&mut model, big, run)
        } else {
            arc_then_arc(&mut model, big, small)
        };
        let pipe = make_pipe_shell_with(
            &mut model,
            &profile,
            &spine,
            &law,
            PipeCorners::Mitre,
            skin,
            T,
        )
        .unwrap_or_else(|e| panic!("{case}: {e}"))
        .shape;
        let found = ogeom::algo::check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{case}: {found}");
        let measured =
            volume_properties(&model, &pipe, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
        assert!(measured.deflection == 0.0, "{case}: not integrated exactly");
        let v = measured.mass;
        let (area, perimeter) = if has_hole { holed } else { round };
        let turned = case == "holed, leaning";
        let want = match (line, has_hole) {
            (true, _) => mitred_against_an_arc(big, run, area, half, disc_width),
            (false, true) => {
                let second: &dyn Fn(f64) -> f64 =
                    if turned { &diamond_width } else { &holed_width };
                mitred_between_arcs(big, small, area, reach, (&holed_width, second))
            }
            (false, false) => {
                mitred_between_arcs(big, small, area, half, (&disc_width, &disc_width))
            }
        };
        let path = if line {
            run
        } else {
            (small + reach) * PI / 2.0
        } + (big + reach) * PI / 2.0
            + 4.0 * reach;
        let bound = skin * perimeter * path;
        assert!(
            (v - want).abs() < bound,
            "{case}: {v} against {want}, bound {bound}"
        );
    }
}

/// A closed spine of two arcs of radius `big` in the `xy` plane, through
/// the origin and `(0, 2 half_chord, 0)`, their centres `(∓d, half_chord,
/// 0)`: a lens, each arc turning anticlockwise, the spine turning sharply
/// at both points. Returned with the angle each arc turns through.
fn lens(model: &mut Model, big: f64, half_chord: f64) -> (Shape, f64) {
    let d = (big * big - half_chord * half_chord).sqrt();
    let spread = half_chord.atan2(d);
    let (a, b) = (
        ogeom::algo::make_vertex(model, Point::ORIGIN).shape,
        ogeom::algo::make_vertex(model, Point::new(0.0, 2.0 * half_chord, 0.0)).shape,
    );
    let mut arc = |x: f64, from: f64, ends: (&Shape, &Shape)| {
        let frame = Frame::new(
            Point::new(x, half_chord, 0.0),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let circle = CircleCurve::new(Circle::new(frame, big, T).unwrap());
        ogeom::algo::make_edge_between(
            model,
            circle.into(),
            (from - spread, from + spread),
            ends.0,
            ends.1,
            T,
        )
        .unwrap()
        .shape
    };
    let first = arc(-d, 0.0, (&a, &b));
    let second = arc(d, PI, (&b, &a));
    (
        make_wire(model, &[first, second], T).unwrap().shape,
        2.0 * spread,
    )
}

/// A disc swept round [`lens`] under the level binormal: two skinned runs
/// mitred at both corners, the second run starting on the first's end
/// circle mirrored, and the first run's start meeting the second's end at
/// the wrap. On a planar spine the binormal's frames are level, and the
/// disc is symmetric, so at each corner each run is run on by the wedge
/// `u tan(φ)` at `u` outside the corner and cut by
/// [`cut_round_an_arc`] at `ρ = big - u` inside it, `2φ` the turn there,
/// both arcs turning with the corner. Pappus gives `π` times the length
/// for the runs; the four corners' columns are integrated by the
/// midpoint rule. The volume is integrated on the faces exactly, so it
/// strays by at most the skin's tolerance over the walls: the perimeter
/// times the longest path a point runs, round both arcs at the outer
/// radius and past each corner by its wedge.
#[test]
fn a_disc_round_a_lens_mitres_its_skinned_runs_at_both_corners() {
    const COLUMNS: u32 = 100_000;
    let (big, half_chord, half, skin) = (15.0, 10.0, 1.0, 1e-3);
    let mut model = Model::new();
    let (spine, turned) = lens(&mut model, big, half_chord);
    // The first arc leaves the origin heading `(sin s, cos s)` for `s` half
    // its turn and comes into the far point heading `(-sin s, cos s)`, where
    // the second leaves heading `(-sin s, -cos s)`: each corner turns
    // through `π - turned`.
    let half_turn = (PI - turned) / 2.0;
    let start = Vector::new((turned / 2.0).sin(), (turned / 2.0).cos(), 0.0);
    let frame = Frame::new(
        Point::ORIGIN,
        Direction::new(start, T).unwrap(),
        Direction::Z,
        T,
    )
    .unwrap();
    let wire = circle_wire(&mut model, frame, half);
    let profile = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let pipe = make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::Binormal(Direction::Z),
        PipeCorners::Mitre,
        skin,
        T,
    )
    .unwrap()
    .shape;
    let found = ogeom::algo::check(&model, &pipe, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let measured =
        volume_properties(&model, &pipe, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
    assert!(measured.deflection == 0.0, "not integrated exactly");
    let du = half / f64::from(COLUMNS);
    let corner: f64 = (0..COLUMNS)
        .map(|k| {
            let u = (f64::from(k) + 0.5) * du;
            let width = 2.0 * (half * half - u * u).max(0.0).sqrt();
            (u * half_turn.tan() - cut_round_an_arc(big, big - u, half_turn, false)) * width * du
        })
        .sum();
    let want = PI * half * half * 2.0 * big * turned + 4.0 * corner;
    let path = 2.0 * (big + half) * turned + 4.0 * half * half_turn.tan();
    let bound = skin * 2.0 * PI * half * path;
    let v = measured.mass;
    assert!(
        (v - want).abs() < bound,
        "{v} against {want}, bound {bound}"
    );
}

/// The square of side `2 half` about the origin on the plane of `frame`,
/// its sides along the frame's `x` and `y`, with a circular hole of
/// `radius` about the origin where one is given.
fn lens_section(model: &mut Model, frame: Frame, half: f64, radius: Option<f64>) -> Shape {
    let (across, up) = (frame.y().vector(), frame.x().vector());
    let pts = [(-half, -half), (half, -half), (half, half), (-half, half)]
        .map(|(s, t)| Point::ORIGIN + across * s + up * t);
    let mut rings = vec![make_polygon(model, &pts, true, T).unwrap().shape];
    if let Some(radius) = radius {
        rings.push(circle_wire(model, frame, radius));
    }
    make_face(
        model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &rings,
        T,
    )
    .unwrap()
    .shape
}

/// A holed square round [`lens`] under the level and the leaning
/// binormal: two skinned runs mitred at both corners, whose hole walls
/// are the tori of one radius the two arcs carry the hole's circle round,
/// mirror images across each mitre plane. Each run sweeps the hole's
/// circle, centred on the spine and square to it, as the torus it is
/// under any frame, so the holed pipe is the square's pipe under the same
/// law less that tube mitred round the lens. The tube is the disc lens's
/// closed form at the hole's radius: Pappus round the arcs, and at each
/// corner each run's wedge `u tan(φ)` outside it less [`cut_round_an_arc`]
/// at `ρ = big - u` inside it, integrated by the midpoint rule.
#[test]
fn a_holed_square_round_a_lens_mitres_its_tube_under_the_binormal() {
    const COLUMNS: u32 = 100_000;
    let (big, half_chord, half, radius, skin) = (15.0, 10.0, 1.0, 0.4, 1e-3);
    let lean = Direction::new(Vector::new(1.0, 0.0, 1.0), T).unwrap();
    let du = radius / f64::from(COLUMNS);
    for (case, law) in [
        ("level", PipeLaw::Binormal(Direction::Z)),
        ("leaning", PipeLaw::Binormal(lean)),
    ] {
        let sweep = |hole: Option<f64>| {
            let mut model = Model::new();
            let (spine, turned) = lens(&mut model, big, half_chord);
            let start = Vector::new((turned / 2.0).sin(), (turned / 2.0).cos(), 0.0);
            let frame = Frame::new(
                Point::ORIGIN,
                Direction::new(start, T).unwrap(),
                Direction::Z,
                T,
            )
            .unwrap();
            let profile = lens_section(&mut model, frame, half, hole);
            let pipe = make_pipe_shell_with(
                &mut model,
                &profile,
                &spine,
                &law,
                PipeCorners::Mitre,
                skin,
                T,
            )
            .unwrap_or_else(|e| panic!("{case}: {e}"))
            .shape;
            let found = ogeom::algo::check(&model, &pipe, T).unwrap();
            assert!(found.is_valid(), "{case}: {found}");
            let measured =
                volume_properties(&model, &pipe, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
            assert!(measured.deflection == 0.0, "{case}: not integrated exactly");
            (measured.mass, turned)
        };
        let (holed, turned) = sweep(Some(radius));
        let (solid, _) = sweep(None);
        let half_turn = (PI - turned) / 2.0;
        let corner: f64 = (0..COLUMNS)
            .map(|k| {
                let u = (f64::from(k) + 0.5) * du;
                let width = 2.0 * (radius * radius - u * u).max(0.0).sqrt();
                (u * half_turn.tan() - cut_round_an_arc(big, big - u, half_turn, false))
                    * width
                    * du
            })
            .sum();
        let tube = PI * radius * radius * 2.0 * big * turned + 4.0 * corner;
        assert!(
            (solid - holed - tube).abs() < tube * 1e-9,
            "{case}: {solid} less {holed} against the tube's {tube}"
        );
    }
}

/// A face round [`lens`] under `law`, the face standing in the plane of
/// the frame the lens starts with, its `x` up: the square of side `2 half`
/// about the spine where `half` is given, a disc of `radius` standing
/// `off` (across, up) from the spine otherwise, or as a hole in the square
/// where both are given. Swept with `skin`, checked valid and integrated
/// exactly; the volume, and the lens's turn.
fn round_the_lens(
    law: &PipeLaw<'_>,
    half: Option<f64>,
    disc: Option<(f64, (f64, f64))>,
    skin: f64,
) -> (f64, f64) {
    let (big, half_chord) = (15.0, 10.0);
    let mut model = Model::new();
    let (spine, turned) = lens(&mut model, big, half_chord);
    let start = Vector::new((turned / 2.0).sin(), (turned / 2.0).cos(), 0.0);
    let frame = Frame::new(
        Point::ORIGIN,
        Direction::new(start, T).unwrap(),
        Direction::Z,
        T,
    )
    .unwrap();
    let (across, up) = (frame.y().vector(), frame.x().vector());
    let mut rings = Vec::new();
    if let Some(half) = half {
        let pts = [(-half, -half), (half, -half), (half, half), (-half, half)]
            .map(|(s, t)| Point::ORIGIN + across * s + up * t);
        rings.push(make_polygon(&mut model, &pts, true, T).unwrap().shape);
    }
    if let Some((radius, (s, t))) = disc {
        let centre =
            Frame::new(Point::ORIGIN + across * s + up * t, frame.z(), frame.x(), T).unwrap();
        rings.push(circle_wire(&mut model, centre, radius));
    }
    let profile = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &rings,
        T,
    )
    .unwrap()
    .shape;
    let pipe = make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        law,
        PipeCorners::Mitre,
        skin,
        T,
    )
    .unwrap_or_else(|e| panic!("{half:?} {disc:?}: {e}"))
    .shape;
    let found = ogeom::algo::check(&model, &pipe, T).unwrap();
    assert!(found.is_valid(), "{half:?} {disc:?}: {found}");
    let measured =
        volume_properties(&model, &pipe, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
    assert!(
        measured.deflection == 0.0,
        "{half:?} {disc:?}: not integrated exactly"
    );
    (measured.mass, turned)
}

/// Sections off the spine round [`lens`], two skinned runs mitred at both
/// corners.
///
/// - Under the level binormal a planar spine's frames are its
///   rotation-minimizing ones, turned across each corner about the
///   corner's axis as that frame is, so a disc off the spine sweeps the
///   pipe the rotation-minimizing frame does: arcs revolved and mitred,
///   built apart. The two agree within the skin's tolerance over the
///   wall: the perimeter times the longest path a point runs.
/// - A hole raised off the spine, square to the lens's plane, keeps its
///   height under those frames, so a square with it under either law is
///   the square's pipe less the mitred tube about the spine of
///   [`a_holed_square_round_a_lens_mitres_its_tube_under_the_binormal`],
///   within the skin's tolerance over the hole's wall.
/// - Under the leaning binormal a square with a hole off the spine is the
///   square's pipe less the hole's own, swept alone under the same law.
#[test]
fn sections_off_the_spine_round_a_lens_mitre_their_skinned_runs() {
    const COLUMNS: u32 = 100_000;
    let (big, half, radius, skin) = (15.0, 1.0, 0.4, 1e-3);
    let lean = Direction::new(Vector::new(1.0, 0.0, 1.0), T).unwrap();
    let level = PipeLaw::Binormal(Direction::Z);
    let carried = PipeLaw::RotationMinimizing;

    let off = (0.3, 0.4);
    let (binormal, turned) = round_the_lens(&level, None, Some((half, off)), skin);
    let (rotation_minimizing, _) = round_the_lens(&carried, None, Some((half, off)), skin);
    let half_turn = (PI - turned) / 2.0;
    let path = |r: f64| 2.0 * (big + r) * turned + 4.0 * r * half_turn.tan();
    let reach = half + off.0.hypot(off.1);
    let bound = skin * 2.0 * PI * half * path(reach);
    assert!(
        (binormal - rotation_minimizing).abs() < bound,
        "{binormal} against {rotation_minimizing}, bound {bound}"
    );

    let du = radius / f64::from(COLUMNS);
    let corner: f64 = (0..COLUMNS)
        .map(|k| {
            let u = (f64::from(k) + 0.5) * du;
            let width = 2.0 * (radius * radius - u * u).max(0.0).sqrt();
            (u * half_turn.tan() - cut_round_an_arc(big, big - u, half_turn, false)) * width * du
        })
        .sum();
    let tube = PI * radius * radius * 2.0 * big * turned + 4.0 * corner;
    let bound = skin * 2.0 * PI * radius * path(radius);
    for (case, law) in [("level", &level), ("carried", &carried)] {
        let (square, _) = round_the_lens(law, Some(half), None, skin);
        let (holed, _) = round_the_lens(law, Some(half), Some((radius, (0.0, 0.5))), skin);
        assert!(
            (square - holed - tube).abs() < bound,
            "{case}: {square} less {holed} against the tube's {tube}, bound {bound}"
        );
    }

    let law = PipeLaw::Binormal(lean);
    let (square, _) = round_the_lens(&law, Some(half), None, skin);
    let (holed, _) = round_the_lens(&law, Some(half), Some((radius, off)), skin);
    let (hole, _) = round_the_lens(&law, None, Some((radius, off)), skin);
    let bound = skin * 2.0 * PI * radius * path(radius + off.0.hypot(off.1));
    assert!(
        (square - holed - hole).abs() < bound,
        "{square} less {holed} against the hole's {hole}, bound {bound}"
    );
}
