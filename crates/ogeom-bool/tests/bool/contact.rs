//! Faces lying on each other's surfaces, where one's boundary runs along
//! the other's.

#![allow(clippy::unwrap_used)]

use ogeom_core::Tolerances;
use ogeom_geom::{Curve, LineCurve};
use ogeom_math::{Frame, Point, Vector};

const T: Tolerances = Tolerances::millimetres();

/// A cube two on a side, and a prism one high standing on its floor
/// against its wall `y = 0`, its foot running along that wall from
/// `x = 0.5` to `x = 1.9` and then a further `short` as a fitted edge:
/// the cube less the prism keeps eight less the prism's area. Measured
/// along the cube's edge at its own spacing, a fitted edge along it
/// shorter than three spacings was never along it, and its stretch of
/// the cube's edge split the prism's floor again as a section running
/// back along the floor's own boundary: the floor kept no piece.
#[test]
fn a_short_fitted_edge_along_the_other_solid_s_edge_is_boundary() {
    for short in [0.0005, 0.005, 0.05] {
        let mut model = ogeom_topo::Model::new();
        let cube = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T)
            .unwrap()
            .shape;
        let foot = 1.9 + short;
        let corners = [
            Point::new(0.5, 0.0, 0.0),
            Point::new(1.9, 0.0, 0.0),
            Point::new(foot, 0.0, 0.0),
            Point::new(foot, 1.0, 0.0),
            Point::new(0.5, 1.0, 0.0),
        ];
        let vertices: Vec<_> = corners
            .iter()
            .map(|&p| ogeom_algo::make_vertex(&mut model, p).shape)
            .collect();
        let mut edges = Vec::new();
        for k in 0..corners.len() {
            let next = (k + 1) % corners.len();
            let (a, b) = (corners[k], corners[next]);
            let line = Curve::Line(LineCurve::segment(a, b, T).unwrap());
            let (curve, range) = if k == 1 {
                let fitted = line.to_bspline_over((0.0, a.distance(b)), T).unwrap();
                (Curve::BSpline(fitted), (0.0, 1.0))
            } else {
                (line, (0.0, a.distance(b)))
            };
            edges.push(
                ogeom_algo::make_edge_between(
                    &mut model,
                    curve,
                    range,
                    &vertices[k],
                    &vertices[next],
                    T,
                )
                .unwrap()
                .shape,
            );
        }
        let wire = ogeom_algo::make_wire(&mut model, &edges, T).unwrap().shape;
        let plane = ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(Frame::WORLD));
        let floor = ogeom_algo::make_face(&mut model, plane.into(), &[wire], T)
            .unwrap()
            .shape;
        let prism = ogeom_algo::make_prism(&mut model, &floor, Vector::new(0.0, 0.0, 1.0), T)
            .unwrap()
            .shape;
        let cut = ogeom_bool::cut(&mut model, &cube, &prism, T).unwrap().shape;
        assert!(ogeom_algo::check(&model, &cut, T).unwrap().is_valid());
        let fine = ogeom_mesh::Deflection {
            chord: 1e-4,
            ..ogeom_mesh::Deflection::default()
        };
        let volume = ogeom_algo::volume_properties(&model, &cut, fine, T)
            .unwrap()
            .mass;
        let want = 8.0 - (foot - 0.5);
        assert!(
            (volume - want).abs() < 1e-9,
            "{short}: {volume} against {want}"
        );
    }
}

/// A pad on the unit circle and a pad on an ellipse one semi-axis `1 + w`
/// long, the other `1`, stacked face to face: the ellipse's cap holds the
/// circle's and passes beyond it in two lunes, each tangent to the circle
/// at both ends and at most `w` wide, far thinner than the band the
/// boolean reads the other solid within. The fuse keeps both pads whole,
/// `π (2 + w)`. Every probe of a lune read it on the circle's cap, and
/// the fuse refused.
#[test]
fn lunes_beside_a_cap_tangent_to_them_are_kept() {
    use ogeom_geom::{CircleCurve, Curve3d as _, EllipseCurve, PlaneSurface};
    use ogeom_math::{Circle, Ellipse, Plane};
    let pad = |model: &mut ogeom_topo::Model, curve: Curve, height: f64| {
        let range = curve.domain();
        let edge = ogeom_algo::make_edge(model, curve, range, T).unwrap().shape;
        let wire = ogeom_algo::make_wire(model, &[edge], T).unwrap().shape;
        let plane = PlaneSurface::over(Plane::XY, (-5.0, 5.0), (-5.0, 5.0)).unwrap();
        let face = ogeom_algo::make_face(model, plane.into(), &[wire], T)
            .unwrap()
            .shape;
        ogeom_algo::make_prism(model, &face, Vector::new(0.0, 0.0, height), T)
            .unwrap()
            .shape
    };
    for w in [2e-5, 3e-5] {
        for circle_first in [true, false] {
            let mut model = ogeom_topo::Model::new();
            let circle =
                Curve::Circle(CircleCurve::new(Circle::new(Frame::WORLD, 1.0, T).unwrap()));
            let ellipse = Curve::Ellipse(EllipseCurve::new(
                Ellipse::new(Frame::WORLD, 1.0 + w, 1.0, T).unwrap(),
            ));
            let below = pad(&mut model, circle, -1.0);
            let above = pad(&mut model, ellipse, 1.0);
            let (a, b) = if circle_first {
                (&below, &above)
            } else {
                (&above, &below)
            };
            let fused = ogeom_bool::fuse(&mut model, a, b, T).unwrap().shape;
            assert!(ogeom_algo::check(&model, &fused, T).unwrap().is_valid());
            let fine = ogeom_mesh::Deflection {
                chord: 1e-4,
                ..ogeom_mesh::Deflection::default()
            };
            let volume = ogeom_algo::volume_properties(&model, &fused, fine, T)
                .unwrap()
                .mass;
            let want = core::f64::consts::PI * (2.0 + w);
            assert!(
                (volume - want).abs() < 1e-9,
                "w {w}, circle first {circle_first}: {volume} against {want}"
            );
        }
    }
}

/// A 20 x 20 x 10 block with its edge `y = 0` at the top rounded at one by
/// cutting the wedge between the edge and its rounding cylinder, then cut
/// by the same wedge along its edge `x = 20`. Each rounding's rail lies on
/// the other wedge's end face and crosses the circle of that face's arc
/// twice: at the arc's end, and a diameter away, at `(18, 0, 9)` and
/// `(20, 2, 9)`, outside every face of the wedge. Only the first is on the
/// face's boundary. Paved at both, each rail kept a vertex no face needs.
/// The block keeps 4000 less the two wedges plus their overlap at the
/// corner, `5/3 - pi/2`.
#[test]
fn a_contact_splits_only_where_it_crosses_the_target_edge_itself() {
    use core::f64::consts::{FRAC_PI_2, FRAC_PI_4};
    use ogeom_math::Direction;
    use ogeom_topo::{ShapeType, explore_unique};
    let mut model = ogeom_topo::Model::new();
    // The bar from `corner` sized `size`, less the cylinder of radius one
    // along `axis` from `start`, run a unit past the bar at either end so
    // no cap lies on the bar's ends.
    let wedge = |model: &mut ogeom_topo::Model,
                 corner: Point,
                 size: (f64, f64, f64),
                 start: Point,
                 axis: Direction| {
        let at = Frame::new(corner, Direction::Z, Direction::X, T).unwrap();
        let bar = ogeom_algo::make_box(model, at, size, T).unwrap().shape;
        let up = Frame::new(start, axis, Direction::Z.cross(axis, T).unwrap(), T).unwrap();
        let drum = ogeom_algo::make_cylinder(model, up, 1.0, 22.0, T)
            .unwrap()
            .shape;
        ogeom_bool::cut(model, &bar, &drum, T).unwrap().shape
    };
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let along_x = wedge(
        &mut model,
        Point::new(0.0, 0.0, 9.0),
        (20.0, 1.0, 1.0),
        Point::new(-1.0, 1.0, 9.0),
        Direction::X,
    );
    let along_y = wedge(
        &mut model,
        Point::new(19.0, 0.0, 9.0),
        (1.0, 20.0, 1.0),
        Point::new(19.0, -1.0, 9.0),
        Direction::Y,
    );
    let rounded = ogeom_bool::cut(&mut model, &block, &along_x, T)
        .unwrap()
        .shape;
    let both = ogeom_bool::cut(&mut model, &rounded, &along_y, T)
        .unwrap()
        .shape;
    assert!(ogeom_algo::check(&model, &both, T).unwrap().is_valid());
    let volume = ogeom_algo::volume_properties(&model, &both, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass;
    let wedge_volume = 20.0 * (1.0 - FRAC_PI_4);
    let want = 4000.0 - 2.0 * wedge_volume + (5.0 / 3.0 - FRAC_PI_2);
    // The faces are planes and cylinders, measured exactly.
    assert!((volume - want).abs() < 1e-9, "{volume} against {want}");
    for v in explore_unique(&model, &both, ShapeType::Vertex).unwrap() {
        let at = model.node(&v).unwrap().data().as_vertex().unwrap().point;
        for off in [Point::new(18.0, 0.0, 9.0), Point::new(20.0, 2.0, 9.0)] {
            assert!(at.distance(off) > 1e-3, "a vertex off both faces at {at:?}");
        }
    }
}

/// A block `x` in `[-7, 0]` whose wall `y - z/10 = 0.8` leans out over its
/// floor, fused with a drum of radius three two high on the origin, the
/// block's top crease rounded at one half by cutting its spandrel (the
/// section between the top, the wall and the ball's arc, run along `x`
/// from `-8` to `-1`) less the drum. The round's cylinder is tangent to
/// the wall along its contact, and where that line meets the drum the
/// spandrel's edge on the drum (a fitted section of the two cylinders)
/// leaves the drum's own edge along the wall (an ellipse) tangentially.
/// Read as running along the ellipse while the two stood within the
/// fitted width of each other, the fitted edge split the ellipse where
/// they parted by it, and the spandrel's own copy of the ellipse did not
/// split there: the kept pieces did not close. The block keeps all but
/// the spandrel from its end to the drum, the section's integral of
/// `7 - sqrt(9 - y^2)`.
#[test]
fn a_fitted_edge_leaving_an_edge_tangentially_from_its_own_end_runs_beside_it() {
    use ogeom_geom::{CircleCurve, PlaneSurface};
    use ogeom_math::{Circle, Direction, Plane};
    let mut model = ogeom_topo::Model::new();
    let face_through =
        |model: &mut ogeom_topo::Model, corners: &[Point], arc: Option<(usize, Point)>| {
            let vertices: Vec<_> = corners
                .iter()
                .map(|&p| ogeom_algo::make_vertex(model, p).shape)
                .collect();
            let mut edges = Vec::new();
            for k in 0..corners.len() {
                let next = (k + 1) % corners.len();
                let (a, b) = (corners[k], corners[next]);
                let (curve, range) = match arc {
                    Some((at, centre)) if at == k => {
                        let (x, y) = (a - centre, b - centre);
                        let frame = Frame::new(
                            centre,
                            Direction::new(x.cross(y), T).unwrap(),
                            Direction::new(x, T).unwrap(),
                            T,
                        )
                        .unwrap();
                        let circle = Circle::new(frame, x.magnitude(), T).unwrap();
                        (
                            Curve::Circle(CircleCurve::new(circle)),
                            (0.0, x.angle(y, T).unwrap()),
                        )
                    }
                    _ => (
                        Curve::Line(LineCurve::segment(a, b, T).unwrap()),
                        (0.0, a.distance(b)),
                    ),
                };
                edges.push(
                    ogeom_algo::make_edge_between(
                        model,
                        curve,
                        range,
                        &vertices[k],
                        &vertices[next],
                        T,
                    )
                    .unwrap()
                    .shape,
                );
            }
            let wire = ogeom_algo::make_wire(model, &edges, T).unwrap().shape;
            let plane = Plane::through(corners[0], Direction::X);
            ogeom_algo::make_face(model, PlaneSurface::new(plane).into(), &[wire], T)
                .unwrap()
                .shape
        };
    let at = |x: f64, y: f64, z: f64| Point::new(x, y, z);
    let end = face_through(
        &mut model,
        &[
            at(-7.0, -5.0, 0.0),
            at(-7.0, 0.8, 0.0),
            at(-7.0, 1.0, 2.0),
            at(-7.0, -5.0, 2.0),
        ],
        None,
    );
    let block = ogeom_algo::make_prism(&mut model, &end, Vector::new(7.0, 0.0, 0.0), T)
        .unwrap()
        .shape;
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 3.0, 2.0, T)
        .unwrap()
        .shape;
    let solid = ogeom_bool::fuse(&mut model, &block, &drum, T)
        .unwrap()
        .shape;
    // The ball of radius one half in the crease: a half below the top, a
    // half in from the wall along its normal.
    let r = 0.5;
    let lean = 1.01_f64.sqrt();
    let yc = 0.95 - r * lean;
    let (ny, nz) = (1.0 / lean, -0.1 / lean);
    let on_wall = (yc + r * ny, 1.5 + r * nz);
    let section = face_through(
        &mut model,
        &[
            at(-8.0, 1.0, 2.0),
            at(-8.0, yc, 2.0),
            at(-8.0, on_wall.0, on_wall.1),
        ],
        Some((1, at(-8.0, yc, 1.5))),
    );
    let spandrel = ogeom_algo::make_prism(&mut model, &section, Vector::new(7.0, 0.0, 0.0), T)
        .unwrap()
        .shape;
    let below = Frame::new(at(0.0, 0.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let tall = ogeom_algo::make_cylinder(&mut model, below, 3.0, 4.0, T)
        .unwrap()
        .shape;
    let past = ogeom_bool::cut(&mut model, &spandrel, &tall, T)
        .unwrap()
        .shape;
    // Kept below `z = cap`, a little above the contact: the fitted edge
    // short enough that its stretch near the ellipse spans samples.
    let cap = 1.52;
    let slab = Frame::new(at(-9.0, -1.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let slab = ogeom_algo::make_box(&mut model, slab, (10.0, 3.0, cap), T)
        .unwrap()
        .shape;
    let tool = ogeom_bool::common(&mut model, &past, &slab, T)
        .unwrap()
        .shape;
    let rounded = ogeom_bool::cut(&mut model, &solid, &tool, T).unwrap().shape;
    assert!(ogeom_algo::check(&model, &rounded, T).unwrap().is_valid());
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let volume = |model: &ogeom_topo::Model, s: &ogeom_topo::Shape| {
        ogeom_algo::volume_properties(model, s, fine, T)
            .unwrap()
            .mass
    };
    // Across the section at height `z = 1.5 + r sin t`, from the arc to the
    // wall; the length's antiderivative in `y` is closed, the angle taken
    // by Simpson's rule.
    let length = |y: f64| 7.0 * y - (y * (9.0 - y * y).sqrt() / 2.0 + 4.5 * (y / 3.0).asin());
    let across = |t: f64| {
        let z = r.mul_add(t.sin(), 1.5);
        let (from, to) = (r.mul_add(t.cos(), yc), 0.1_f64.mul_add(z, 0.8));
        (length(to) - length(from)) * r * t.cos()
    };
    let (t0, t1) = (nz.atan2(ny), ((cap - 1.5) / r).asin());
    let n = 4000;
    let h = (t1 - t0) / f64::from(n);
    let mut sum = across(t0) + across(t1);
    for k in 1..n {
        sum += across(f64::from(k).mul_add(h, t0)) * if k % 2 == 1 { 4.0 } else { 2.0 };
    }
    let removed = sum * h / 3.0;
    let got = volume(&model, &solid) - volume(&model, &rounded);
    assert!((got - removed).abs() < 1e-7, "{got} against {removed}");
}
