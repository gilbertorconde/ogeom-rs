//! Solids offset by `offset_shape`, their topology kept one for one, held
//! to the volumes their grown sizes give.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

fn torus_volume(major: f64, minor: f64) -> f64 {
    2.0 * PI * PI * major * minor * minor
}

/// A whole torus is one face bounded by its two seams, each walked twice,
/// crossing at one vertex: no ends for the faces' wires to be assembled
/// from, and the offset solid did not close. It is the moved torus whole,
/// out and in, on an oblique axis as on the world's.
#[test]
fn a_whole_torus_offsets_out_and_in() {
    let tilted = Frame::new(
        Point::new(1.0, -2.0, 0.5),
        Direction::new(Vector::new(0.3, -0.5, 0.8), T).unwrap(),
        Direction::new(Vector::new(0.0, 0.8, 0.5), T).unwrap(),
        T,
    )
    .unwrap();
    for frame in [Frame::WORLD, tilted] {
        for d in [0.4, -0.3] {
            let mut model = Model::new();
            let ring = ogeom::algo::make_torus(&mut model, frame, 3.0, 1.0, T)
                .unwrap()
                .shape;
            let grown = ogeom::offset::offset_shape(&mut model, &ring, d, T)
                .unwrap()
                .shape;
            let diagnosis = ogeom::algo::check(&model, &grown, T).unwrap();
            assert!(diagnosis.is_valid(), "offset {d}: {diagnosis}");
            let (got, want) = (volume(&model, &grown), torus_volume(3.0, 1.0 + d));
            assert!(
                (got - want).abs() <= want * 1e-9,
                "offset {d}: {got} against {want}"
            );
        }
    }
}

/// A block with a torus-shaped void: grown, the block's faces move out and
/// the void's tube narrows by as much.
#[test]
fn a_torus_void_narrows_as_its_block_grows() {
    let mut model = Model::new();
    let corner = Frame::new(Point::new(-6.0, -6.0, -3.0), Direction::Z, Direction::X, T).unwrap();
    let block = ogeom::algo::make_box(&mut model, corner, (12.0, 12.0, 6.0), T)
        .unwrap()
        .shape;
    let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, 3.0, 1.0, T)
        .unwrap()
        .shape;
    let holed = ogeom::boolean::cut(&mut model, &block, &ring, T)
        .unwrap()
        .shape;
    let d = 0.25;
    let grown = ogeom::offset::offset_shape(&mut model, &holed, d, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &grown, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let want = (12.0 + 2.0 * d) * (12.0 + 2.0 * d) * (6.0 + 2.0 * d) - torus_volume(3.0, 1.0 - d);
    let got = volume(&model, &grown);
    assert!((got - want).abs() <= want * 1e-9, "{got} against {want}");
}

/// A ball of `radius` about the world's z axis, one face whose seam is two
/// edges leaving the equator vertex on its x axis, one up to the north
/// pole along its meridian's parameter and one down to the south pole
/// against it.
fn ball_seamed_both_ways_from_the_equator(model: &mut Model, radius: f64) -> Shape {
    use ogeom::geom::{CircleCurve, Line2d, PlanarCurve, SphereSurface};
    use ogeom::math::{Axis2, Circle, Direction2, Point2, Sphere, Vector2};
    use ogeom::topo::{EdgeData, Location};
    let tau = 2.0 * PI;
    let half = PI / 2.0;
    let equator = ogeom::algo::make_vertex(model, Point::new(radius, 0.0, 0.0)).shape;
    let north = ogeom::algo::make_vertex(model, Point::new(0.0, 0.0, radius)).shape;
    let south = ogeom::algo::make_vertex(model, Point::new(0.0, 0.0, -radius)).shape;
    // Each meridian's angle is the latitude, upward; turned over, downward.
    let meridian = |normal: Direction| {
        Circle::new(
            Frame::new(Point::ORIGIN, normal, Direction::X, T).unwrap(),
            radius,
            T,
        )
        .unwrap()
    };
    let up = ogeom::algo::make_edge_between(
        model,
        CircleCurve::new(meridian(-Direction::Y)).into(),
        (0.0, half),
        &equator,
        &north,
        T,
    )
    .unwrap()
    .shape;
    let down = ogeom::algo::make_edge_between(
        model,
        CircleCurve::new(meridian(Direction::Y)).into(),
        (0.0, half),
        &equator,
        &south,
        T,
    )
    .unwrap()
    .shape;
    let pole = |model: &mut Model, at: &Shape| {
        let mut data = EdgeData::new();
        data.degenerate = true;
        model.add_edge(data, &[at.clone(), at.clone()]).unwrap()
    };
    let (top, bottom) = (pole(model, &north), pole(model, &south));
    let surface = model
        .geometry_mut()
        .add_surface(SphereSurface::new(Sphere::new(Frame::WORLD, radius, T).unwrap()).into());
    let line = |from: (f64, f64), along: (f64, f64)| -> PlanarCurve {
        let direction = Direction2::new(Vector2::new(along.0, along.1), T).unwrap();
        Line2d::over(
            Axis2::new(Point2::new(from.0, from.1), direction),
            -1.0,
            tau + 1.0,
        )
        .unwrap()
        .into()
    };
    let here = Location::identity();
    for (row, edge) in [(half, &top), (-half, &bottom)] {
        ogeom::algo::attach_pcurve(
            model,
            edge,
            line((0.0, row), (1.0, 0.0)),
            surface,
            here.clone(),
            (0.0, tau),
        )
        .unwrap();
    }
    // Up the far column of the chart and down the near one.
    ogeom::algo::attach_seam(
        model,
        &up,
        line((tau, 0.0), (0.0, 1.0)),
        line((0.0, 0.0), (0.0, 1.0)),
        surface,
        here.clone(),
        (0.0, half),
    )
    .unwrap();
    ogeom::algo::attach_seam(
        model,
        &down,
        line((0.0, 0.0), (0.0, -1.0)),
        line((tau, 0.0), (0.0, -1.0)),
        surface,
        here,
        (0.0, half),
    )
    .unwrap();
    let ring = [
        bottom,
        down.reversed(),
        up.clone(),
        top.reversed(),
        up.reversed(),
        down,
    ];
    let wire = ogeom::algo::make_wire(model, &ring, T).unwrap().shape;
    let face = ogeom::algo::make_face_on(model, surface, &[wire], T)
        .unwrap()
        .shape;
    let shell = ogeom::algo::make_shell(model, &[face]).unwrap().shape;
    ogeom::algo::make_solid(model, &[shell]).unwrap().shape
}

/// Offset out and in, the ball seamed both ways from its equator is the
/// ball of the grown radius: its south seam, running against its
/// meridian's parameter, is rebuilt down to the pole, not three quarters
/// of a turn up and round through rows the sphere's chart does not have.
/// The offset meshes closed, and an octant cut from it leaves seven eighths.
#[test]
fn a_ball_seamed_both_ways_from_its_equator_offsets_closed() {
    let mut model = Model::new();
    let ball = ball_seamed_both_ways_from_the_equator(&mut model, 5.0);
    let diagnosis = ogeom::algo::check(&model, &ball, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let ball_volume = |r: f64| 4.0 / 3.0 * PI * r.powi(3);
    assert!((volume(&model, &ball) - ball_volume(5.0)).abs() <= ball_volume(5.0) * 1e-6);
    for d in [0.5, -0.5] {
        let grown = ogeom::offset::offset_shape(&mut model, &ball, d, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom::algo::check(&model, &grown, T).unwrap();
        assert!(diagnosis.is_valid(), "offset {d}: {diagnosis}");
        let want = ball_volume(5.0 + d);
        let got = volume(&model, &grown);
        assert!(
            (got - want).abs() <= want * 1e-6,
            "offset {d}: {got} against {want}"
        );
        assert!(
            ogeom::mesh::triangulate(&model, &grown, Deflection::default(), T)
                .unwrap()
                .is_closed(),
            "offset {d}: the mesh is open"
        );
        let octant = ogeom::algo::make_box(&mut model, Frame::WORLD, (6.0, 6.0, 6.0), T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &grown, &octant, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom::algo::check(&model, &cut, T).unwrap();
        assert!(diagnosis.is_valid(), "offset {d} cut: {diagnosis}");
        let (got, want) = (volume(&model, &cut), want * 7.0 / 8.0);
        assert!(
            (got - want).abs() <= want * 1e-6,
            "offset {d} cut: {got} against {want}"
        );
        assert!(
            ogeom::mesh::triangulate(&model, &cut, Deflection::default(), T)
                .unwrap()
                .is_closed(),
            "offset {d} cut: the mesh is open"
        );
    }
}
