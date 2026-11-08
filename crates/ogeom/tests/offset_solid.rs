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

/// The ellipsoid with semi-axes 5, 3 and 2: a unit ball scaled, one
/// rational B-spline face bounded by its seam and two poles.
fn ellipsoid(model: &mut Model) -> Shape {
    use ogeom::math::{GeneralTransform, Matrix3};
    let unit = ogeom::algo::make_sphere(model, Frame::WORLD, 1.0, T)
        .unwrap()
        .shape;
    let scale = GeneralTransform {
        linear: Matrix3::from_columns(
            Vector::new(5.0, 0.0, 0.0),
            Vector::new(0.0, 3.0, 0.0),
            Vector::new(0.0, 0.0, 2.0),
        ),
        translation: Vector::ZERO,
    };
    ogeom::algo::general_transformed_shape(model, &unit, &scale, T)
        .unwrap()
        .shape
}

/// The volume inside the parallel surface at `d` of the ellipsoid with
/// semi-axes `(a, b, c)`: a third of the flux of the position through it,
/// `P = X + d n`, integrated over latitude by Simpson's rule and round the
/// axis by the trapezoid rule (exact for a periodic integrand to rounding).
fn ellipsoid_parallel_volume((a, b, c): (f64, f64, f64), d: f64) -> f64 {
    let at = |t: f64, f: f64| {
        let x = Vector::new(a * t.cos() * f.cos(), b * t.cos() * f.sin(), c * t.sin());
        let n = Vector::new(x.x / (a * a), x.y / (b * b), x.z / (c * c));
        x + n * (d / n.magnitude())
    };
    let (nt, nf) = (2000, 256);
    let h = 1e-6;
    let mut total = 0.0;
    for i in 0..=nt {
        let t = -PI / 2.0 + PI * f64::from(i) / f64::from(nt);
        let w = if i == 0 || i == nt {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let mut ring = 0.0;
        for j in 0..nf {
            let f = 2.0 * PI * f64::from(j) / f64::from(nf);
            let p = at(t, f);
            let pt = (at(t + h, f) - at(t - h, f)) / (2.0 * h);
            let pf = (at(t, f + h) - at(t, f - h)) / (2.0 * h);
            ring += p.dot(pf.cross(pt));
        }
        total += w * ring;
    }
    total * (PI / f64::from(nt) / 3.0) * (2.0 * PI / f64::from(nf)) / 3.0
}

/// The widest a solid's faces state themselves, and their total area.
fn stated_breadth(model: &Model, shape: &Shape) -> (f64, f64) {
    use ogeom::topo::{ShapeType, explore_unique};
    let mut widest: f64 = 0.0;
    let mut area = 0.0;
    for face in explore_unique(model, shape, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        widest = widest.max(data.tolerance.get());
        area +=
            ogeom::algo::surface_properties(model, &face, Deflection::with_chord(1e-3).unwrap(), T)
                .unwrap()
                .mass;
    }
    (widest, area)
}

/// An ellipsoid, a rational B-spline face with a seam and two poles,
/// offsets out and in: its face moves to a fit of its parallel, its poles
/// along the limit of the normals there. The volume is the parallel
/// surface's own, integrated in closed form, within the fit's stated
/// breadth over the area; the bounds are the semi-axes moved.
#[test]
fn an_ellipsoid_offsets_out_and_in() {
    let axes = (5.0, 3.0, 2.0);
    for d in [1.0, -0.5] {
        let mut model = Model::new();
        let solid = ellipsoid(&mut model);
        let moved = ogeom::offset::offset_shape(&mut model, &solid, d, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom::algo::check(&model, &moved, T).unwrap();
        assert!(diagnosis.is_valid(), "offset {d}: {diagnosis}");
        let (breadth, area) = stated_breadth(&model, &moved);
        assert!(breadth < 1e-3, "offset {d}: the fit states {breadth}");
        let want = ellipsoid_parallel_volume(axes, d);
        let got = volume(&model, &moved);
        eprintln!("offset {d}: volume {got} against {want}, breadth {breadth} over {area}");
        assert!(
            (got - want).abs() <= breadth * area + want * 1e-9,
            "offset {d}: {got} against {want}, breadth {breadth} over {area}"
        );
        let bounds = ogeom::algo::tight_bounds(&model, &moved, T).unwrap();
        let (lo, hi) = (bounds.low().unwrap(), bounds.high().unwrap());
        for (k, half) in [axes.0, axes.1, axes.2].into_iter().enumerate() {
            let reach = half + d;
            let (l, h) = [(lo.x, hi.x), (lo.y, hi.y), (lo.z, hi.z)][k];
            assert!(
                (l + reach).abs() <= 1e-3 && (h - reach).abs() <= 1e-3,
                "offset {d}: axis {k} spans [{l}, {h}] against ±{reach}"
            );
        }
    }
}

/// Inward by 1, past the ellipsoid's smallest radius of curvature (`c² /
/// a` = 0.8 at the ends of its long axis), the parallel folds over itself
/// and the offset is refused by name.
#[test]
fn an_ellipsoid_refuses_an_offset_past_its_curvature() {
    let mut model = Model::new();
    let solid = ellipsoid(&mut model);
    let refused = ogeom::offset::offset_shape(&mut model, &solid, -1.0, T)
        .unwrap_err()
        .to_string();
    assert!(refused.contains("radius of curvature"), "{refused}");
}

/// A prism over a face bounded by a cubic B-spline arch and the line
/// under it: its arched side is an extrusion. Grown, that side moves to
/// the arch's parallel, met by the moved floor, roof and flat side past
/// its old ends; shrunk back, the prism returns.
#[test]
fn a_spline_sided_prism_offsets_out_and_back() {
    use ogeom::geom::{BSplineCurve, Curve, Curve3d as _, LineCurve, PlaneSurface};
    use ogeom::math::{KnotVector, Plane};
    let mut model = Model::new();
    let (a, b) = (Point::ORIGIN, Point::new(4.0, 0.0, 0.0));
    let arch = BSplineCurve::new(
        KnotVector::clamped_uniform(3, 4).unwrap(),
        vec![a, Point::new(1.0, 2.0, 0.0), Point::new(3.0, 2.0, 0.0), b],
        T,
    )
    .unwrap();
    let va = ogeom::algo::make_vertex(&mut model, a).shape;
    let vb = ogeom::algo::make_vertex(&mut model, b).shape;
    let top = ogeom::algo::make_edge_between(
        &mut model,
        Curve::BSpline(arch.clone()),
        (0.0, 1.0),
        &va,
        &vb,
        T,
    )
    .unwrap()
    .shape;
    let floor = LineCurve::segment(b, a, T).unwrap();
    let under =
        ogeom::algo::make_edge_between(&mut model, Curve::from(floor), (0.0, 4.0), &vb, &va, T)
            .unwrap()
            .shape;
    let wire = ogeom::algo::make_wire(&mut model, &[top, under], T)
        .unwrap()
        .shape;
    let ground = Plane::new(Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap());
    let plane = PlaneSurface::over(ground, (-10.0, 10.0), (-10.0, 10.0)).unwrap();
    let face = ogeom::algo::make_face(&mut model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let prism = ogeom::algo::make_prism(&mut model, &face, Vector::new(0.0, 0.0, 3.0), T)
        .unwrap()
        .shape;
    let before = volume(&model, &prism);

    let grown = ogeom::offset::offset_shape(&mut model, &prism, 0.5, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &grown, T).unwrap();
    assert!(diagnosis.is_valid(), "grown: {diagnosis}");
    // Halfway up, the arch's parallel lies on the grown solid's boundary.
    let (breadth, _) = stated_breadth(&model, &grown);
    let boundary = ogeom::algo::PointDistance::of(&model, &grown, T).unwrap();
    for k in 0..=32 {
        let t = f64::from(k) / 32.0;
        let p = arch.point_at(t, T).unwrap();
        let along = arch.d1_at(t, T).unwrap();
        let out = Vector::new(-along.y, along.x, 0.0) * (1.0 / along.magnitude());
        let on = Point::new(p.x, p.y, 1.5) + out * 0.5;
        let (gap, _) = boundary.nearest(on, T).unwrap();
        assert!(
            gap <= breadth.max(1e-7),
            "the parallel at {t} stands {gap} off"
        );
    }

    let back = ogeom::offset::offset_shape(&mut model, &grown, -0.5, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &back, T).unwrap();
    assert!(diagnosis.is_valid(), "shrunk back: {diagnosis}");
    let got = volume(&model, &back);
    eprintln!(
        "grown {} (breadth {breadth}), shrunk back {got} against {before}",
        volume(&model, &grown)
    );
    assert!(
        (got - before).abs() <= before * 1e-3,
        "shrunk back to {got} against {before}"
    );
}
