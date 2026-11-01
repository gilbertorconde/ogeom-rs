//! The additive blends: concave edges, and the revolved seats in every sign.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Point};
use ogeom_topo::{Filter, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

fn fine() -> ogeom_mesh::Deflection {
    ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    }
}

fn volume(model: &ogeom_topo::Model, shape: &ogeom_topo::Shape) -> f64 {
    ogeom_algo::volume_properties(model, shape, fine(), T)
        .unwrap()
        .mass
}

/// An L: a box with a notch, and its re-entrant edge along y at (1, 1).
fn l_bracket(model: &mut ogeom_topo::Model) -> (ogeom_topo::Shape, ogeom_topo::Shape) {
    let block = ogeom_algo::make_box(model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let seat = Frame::new(
        Point::new(1.0, -0.5, 1.0),
        ogeom_math::Direction::Z,
        ogeom_math::Direction::X,
        T,
    )
    .unwrap();
    let notch = ogeom_algo::make_box(model, seat, (2.0, 3.0, 2.0), T).unwrap();
    let cut = ogeom_bool::cut(model, &block.shape, &notch.shape, T).unwrap();
    let edge = explore(model, &cut.shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(model, e)
                .unwrap()
                .is_some_and(|(a, b)| {
                    let p = |v: &ogeom_topo::Shape| {
                        model
                            .node(v)
                            .and_then(|n| n.data().as_vertex().map(|d| d.point))
                            .unwrap()
                    };
                    let (pa, pb) = (p(&a), p(&b));
                    (pa.x - 1.0).abs() < 1e-9
                        && (pa.z - 1.0).abs() < 1e-9
                        && (pb.x - 1.0).abs() < 1e-9
                        && (pb.z - 1.0).abs() < 1e-9
                })
        })
        .expect("the L has its re-entrant edge");
    (cut.shape, edge)
}

#[test]
fn a_concave_chamfer_fills_the_corner_with_its_bevel() {
    let mut model = ogeom_topo::Model::new();
    let (bracket, edge) = l_bracket(&mut model);
    let d = 0.25;
    let result = ogeom_fillet::chamfer_edge(&mut model, &bracket, &edge, d, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let exact = 6.0 + d * d / 2.0 * 2.0;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - exact).abs() < 1e-9,
        "concave chamfer volume {measured} against {exact}"
    );
    assert!(result.history.is_deleted(&edge));
}

/// The rim of a hole: the top edge of a tube's bore.
#[test]
fn a_hole_rim_gains_the_mirrored_toroidal_blend() {
    let mut model = ogeom_topo::Model::new();
    let outer = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 2.0, T).unwrap();
    let bore = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 2.0, T).unwrap();
    let tube = ogeom_bool::cut(&mut model, &outer.shape, &bore.shape, T).unwrap();

    // The bore's top rim: radius 1 at z = 2.
    let edge = explore(&model, &tube.shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(&model, e)
                .unwrap()
                .is_some_and(|(a, _)| {
                    model
                        .node(&a)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .is_some_and(|p| {
                            (p.z - 2.0).abs() < 1e-9 && (p.x.hypot(p.y) - 1.0).abs() < 1e-6
                        })
                })
        })
        .expect("the tube has its bore rim");

    let r = 0.3;
    let result = ogeom_fillet::fillet_edge(&mut model, &tube.shape, &edge, r, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // Pappus over the mirrored meridian cusp at the bore radius.
    let pi = core::f64::consts::PI;
    let hole = 1.0;
    let removed = 2.0
        * pi
        * (hole * r * r + r * r * r / 2.0 - (hole + r) * pi * r * r / 4.0 + r * r * r / 3.0);
    let exact = pi * (4.0 - 1.0) * 2.0 - removed;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - exact).abs() < 2e-3,
        "hole rim fillet volume {measured} against {exact}"
    );
    assert!(result.history.is_deleted(&edge));
}

/// The same mouth in a plate, drilled right through: the bore's rim rounded
/// over where the drill was longer than the plate is thick.
///
/// The drill's length is the condition. A bore exactly as deep as the
/// plate leaves its mouth on the drill cylinder's own rim, and that circle
/// starts where the wall's chart does; a deeper drill leaves it as a
/// section, starting wherever the cut split it. The blend's wedge shares
/// the wall with the plate, and the shared rim's image on that chart is a
/// straight pcurve running from a quarter turn round to a quarter turn
/// past the seam. Two of its samples stand three quarters of a turn apart,
/// which the boolean's tear remover must not read as a wrap: turned back
/// on itself, the pcurve never splits the wall, the wall's strip inside
/// the wedge survives, and the shell does not close.
#[test]
fn a_bore_mouth_in_a_plate_rounds_over_however_deep_the_drill() {
    let pi = core::f64::consts::PI;
    let (plate, bore, r) = (10.0_f64, 4.0_f64, 1.0_f64);
    for deeper in [false, true] {
        let mut model = ogeom_topo::Model::new();
        let block = ogeom_algo::make_box(
            &mut model,
            Frame::new(
                Point::new(-plate, -plate, 0.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            (plate * 2.0, plate * 2.0, 10.0),
            T,
        )
        .unwrap();
        let (low, tall) = if deeper { (-1.0, 12.0) } else { (0.0, 10.0) };
        let seat = Frame::new(
            Point::new(0.0, 0.0, low),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let drill = ogeom_algo::make_cylinder(&mut model, seat, bore, tall, T).unwrap();
        let holed = ogeom_bool::cut(&mut model, &block.shape, &drill.shape, T).unwrap();

        // The mouth, asked for away from the wall's seam so the seam's own
        // edge cannot answer instead.
        let mouth = explore(&model, &holed.shape, Filter::OfType(ShapeType::Edge))
            .unwrap()
            .into_iter()
            .find(|e| {
                use ogeom_geom::Curve3d as _;
                let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = model
                    .node(e)
                    .and_then(|n| n.data().as_edge())
                    .and_then(ogeom_topo::EdgeData::curve3d)
                else {
                    return false;
                };
                let Some(c) = model.geometry().curve(*curve) else {
                    return false;
                };
                (0..=8).all(|i| {
                    let t = range.0 + (range.1 - range.0) * f64::from(i) / 8.0;
                    c.point_at(t, T).is_ok_and(|p| {
                        (p.z - 10.0).abs() < 1e-9 && (p.x.hypot(p.y) - bore).abs() < 1e-6
                    })
                })
            })
            .expect("the plate has its bore mouth");

        let built = ogeom_fillet::fillet_edge(&mut model, &holed.shape, &mouth, r, T)
            .unwrap_or_else(|e| panic!("drill deeper {deeper}: {e}"));
        let diagnosis = ogeom_algo::check(&model, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        assert_eq!(
            explore(&model, &built.shape, Filter::OfType(ShapeType::Face))
                .unwrap()
                .len(),
            8,
            "drill deeper {deeper}: six walls, the bore and the band"
        );
        // Pappus over the mirrored meridian cusp at the bore's radius, the
        // same figure the tube's rim answers to.
        let removed = 2.0
            * pi
            * (bore * r * r + r * r * r / 2.0 - (bore + r) * pi * r * r / 4.0 + r * r * r / 3.0);
        let want = plate * plate * 4.0 * 10.0 - pi * bore * bore * 10.0 - removed;
        let got = volume(&model, &built.shape);
        assert!(
            (got - want).abs() < want * 1e-5,
            "drill deeper {deeper}: {got} against {want}"
        );
    }
}

/// The base of a boss: the concave circular seat where a cylinder stands on
/// a plate.
#[test]
fn a_boss_base_gains_an_additive_toroidal_blend() {
    let mut model = ogeom_topo::Model::new();
    let plate = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 2.0), T).unwrap();
    let seat = Frame::new(
        Point::new(5.0, 5.0, 2.0),
        ogeom_math::Direction::Z,
        ogeom_math::Direction::X,
        T,
    )
    .unwrap();
    let boss = ogeom_algo::make_cylinder(&mut model, seat, 1.0, 3.0, T).unwrap();
    let joined = ogeom_bool::fuse(&mut model, &plate.shape, &boss.shape, T).unwrap();

    // The boss's base circle: radius 1 at z = 2 about (5, 5).
    let edge = explore(&model, &joined.shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(&model, e)
                .unwrap()
                .is_some_and(|(a, _)| {
                    model
                        .node(&a)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .is_some_and(|p| {
                            (p.z - 2.0).abs() < 1e-9
                                && ((p.x - 5.0).hypot(p.y - 5.0) - 1.0).abs() < 1e-6
                        })
                })
        })
        .expect("the joined part has the boss base circle");

    let r = 0.3;
    let result = ogeom_fillet::fillet_edge(&mut model, &joined.shape, &edge, r, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let pi = core::f64::consts::PI;
    let boss_r = 1.0;
    let added = 2.0
        * pi
        * (boss_r * r * r + r * r * r / 2.0 - (boss_r + r) * pi * r * r / 4.0 + r * r * r / 3.0);
    let exact = 200.0 + pi * 3.0 + added;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - exact).abs() < 2e-3,
        "boss base fillet volume {measured} against {exact}"
    );
    assert!(result.history.is_deleted(&edge));
}

/// A plate whose top is `z = 0`, fused with a block over it whose side
/// face leans 15 degrees and leaves the top along a concave seam: the
/// block's top edge, the far side of the narrow strip of that face above
/// the plate, rises at `rise` degrees through a point `lift` over the
/// origin. The seam runs from `ends.0` to `ends.1` along it, measured from
/// the origin; the plate is a prism over a quadrilateral whose ends cross
/// the seam there, turned about z by `turns` degrees; `post` adds a small
/// block standing on the plate in front of the strip.
struct Seat {
    rise: f64,
    lift: f64,
    ends: (f64, f64),
    turns: (f64, f64),
    post: bool,
}

/// The slanted seam the rolling-round tests round: the ball, radius 0.5,
/// touches the plate's top and rests on the block's top edge all along.
const SLANT: Seat = Seat {
    rise: 20.0,
    lift: 0.0,
    ends: (0.6, 1.3),
    turns: (-20.0, 5.0),
    post: false,
};

struct Slant {
    shape: ogeom_topo::Shape,
    seam: ogeom_topo::Shape,
    /// The block's top edge: its direction through the origin.
    far: ogeom_math::Vector,
    /// Outward normals of the block's side face and top face.
    side: ogeom_math::Vector,
    top: ogeom_math::Vector,
    /// The seam's unit direction, and its ends' distances from the origin.
    along: ogeom_math::Vector,
    ends: (f64, f64),
    /// A point of the top edge.
    lift: Point,
}

fn slanted_seam(model: &mut ogeom_topo::Model, seat: &Seat) -> Slant {
    use ogeom_math::{Direction, Plane, Vector};
    let ends = seat.ends;
    let lift = Point::new(0.0, 0.0, seat.lift);
    let (rise, lean) = (seat.rise.to_radians(), 15_f64.to_radians());
    let far = Vector::new(rise.cos(), 0.0, rise.sin());
    let (y0, z0) = (Vector::Y, Vector::new(-rise.sin(), 0.0, rise.cos()));
    let depth = y0 * lean.cos() + z0 * lean.sin();
    let up = z0 * lean.cos() - y0 * lean.sin();
    let block = ogeom_algo::make_box(
        model,
        Frame::new(
            lift - up * 3.0 - far * 2.0,
            Direction::new(up, T).unwrap(),
            Direction::new(far, T).unwrap(),
            T,
        )
        .unwrap(),
        (6.0, 3.0, 3.0),
        T,
    )
    .unwrap();
    let side = -depth;
    let along = {
        let a = Vector::Z.cross(side);
        let a = a * (1.0 / a.magnitude());
        if a.dot(far) < 0.0 { -a } else { a }
    };
    // Where the seam crosses the plane x = 0: on the side face's plane
    // through the top edge, at z = 0.
    let start = Point::new(0.0, side.z * seat.lift / side.y, 0.0);
    // The plate's ends turned about z so each end of the seam keeps the
    // round's foot on the plate's top.
    let corner = |at: f64, angle: f64, y: f64| {
        let through = start + along * at;
        let (s, c) = angle.to_radians().sin_cos();
        let d = Vector::new(-s, c, 0.0);
        through + d * ((y - through.y) / d.y)
    };
    let quad = [
        corner(ends.0, seat.turns.0, -2.0),
        corner(ends.1, seat.turns.1, -2.0),
        corner(ends.1, seat.turns.1, 1.0),
        corner(ends.0, seat.turns.0, 1.0),
    ];
    let wire = ogeom_algo::make_polygon(model, &quad, true, T)
        .unwrap()
        .shape;
    let top = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(model, top.into(), &[wire], T)
        .unwrap()
        .shape;
    let plate = ogeom_algo::make_prism(model, &face, Vector::new(0.0, 0.0, -3.0), T)
        .unwrap()
        .shape;
    let mut shape = ogeom_bool::fuse(model, &block.shape, &plate, T)
        .unwrap()
        .shape;
    if seat.post {
        let post = ogeom_algo::make_box(
            model,
            Frame::new(Point::new(0.8, -1.05, -0.1), Direction::Z, Direction::X, T).unwrap(),
            (0.2, 0.15, 0.4),
            T,
        )
        .unwrap();
        shape = ogeom_bool::fuse(model, &shape, &post.shape, T)
            .unwrap()
            .shape;
    }
    let middle = start + along * ((ends.0 + ends.1) / 2.0);
    let seam = explore(model, &shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(model, e)
                .unwrap()
                .is_some_and(|(a, b)| {
                    let p = |v: &ogeom_topo::Shape| {
                        model
                            .node(v)
                            .and_then(|n| n.data().as_vertex().map(|d| d.point))
                            .unwrap()
                    };
                    // The seam's ends stand where the plate's ends cross it,
                    // to rounding.
                    p(&a).midpoint(p(&b)).distance(middle) < 1e-9
                })
        })
        .expect("the fused part has the seam");
    Slant {
        shape,
        seam,
        far,
        side,
        top: up,
        along,
        ends,
        lift,
    }
}

/// Gauss-Legendre nodes and weights on `[-1, 1]`, by Newton's method on
/// the Legendre polynomial.
fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    use core::f64::consts::PI;
    (0..n)
        .map(|i| {
            #[allow(clippy::cast_precision_loss)]
            let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            let mut dp = 0.0;
            for _ in 0..100 {
                let (mut p0, mut p1) = (1.0, x);
                for k in 2..=n {
                    #[allow(clippy::cast_precision_loss)]
                    let k = k as f64;
                    let p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k;
                    p0 = p1;
                    p1 = p2;
                }
                #[allow(clippy::cast_precision_loss)]
                let nf = n as f64;
                dp = nf * (x * p1 - p0) / (x * x - 1.0);
                let step = p1 / dp;
                x -= step;
                if step.abs() < 1e-16 {
                    break;
                }
            }
            (x, 2.0 / ((1.0 - x * x) * dp * dp))
        })
        .collect()
}

/// The ball of `radius` resting on the plate's top and on the block's top
/// edge: its centre at angle `a` round the edge from the side face's
/// normal toward the top face's, a radius off the plate's top, and the
/// centre's first and second derivatives in `a`.
fn rolling_centre(
    slant: &Slant,
    radius: f64,
    a: f64,
) -> (Point, ogeom_math::Vector, ogeom_math::Vector) {
    use ogeom_math::Vector;
    let (u, v) = (slant.side, slant.top);
    let k = slant.far.dot(Vector::Z);
    let (au, av) = (u.z, v.z);
    let (s, c) = a.sin_cos();
    let reach = (u * c + v * s) * radius;
    let foot = (radius - reach.z) / k;
    let d_foot = radius * (s * au - c * av) / k;
    let dd_foot = radius * (c * au + s * av) / k;
    (
        Point::ORIGIN + slant.far * foot + reach,
        slant.far * d_foot + (v * c - u * s) * radius,
        slant.far * dd_foot - reach,
    )
}

/// Where along the seam the plane square to the centre's path at angle
/// `a` crosses it.
fn rolling_crossing(slant: &Slant, radius: f64, a: f64) -> f64 {
    let (c, d, _) = rolling_centre(slant, radius, a);
    d.dot(c - Point::ORIGIN) / d.dot(slant.along)
}

/// The angles at which the ball's sections cross the seam's two ends.
fn rolling_span(slant: &Slant, radius: f64) -> (f64, f64) {
    let at = |target: f64| {
        let (mut lo, mut hi) = (0.0, core::f64::consts::FRAC_PI_2);
        for _ in 0..200 {
            let m = (lo + hi) / 2.0;
            let f = |a: f64| rolling_crossing(slant, radius, a) - target;
            if (f(lo) > 0.0) == (f(m) > 0.0) {
                lo = m;
            } else {
                hi = m;
            }
        }
        (lo + hi) / 2.0
    };
    (at(slant.ends.0), at(slant.ends.1))
}

/// The volume between the seam's corner and the rolling ball, swept by
/// the ball's sections: each section, in the plane square to the centre's
/// path, is the triangle from the seam to the two contacts less the
/// circular segment the ball's arc cuts off it, and a plane turning with
/// the path sweeps `(|c'| - n'.(p - c)) dA da` at each point `p` of it.
fn rolling_fill(slant: &Slant, radius: f64) -> f64 {
    let (a0, a1) = rolling_span(slant, radius);
    let (lo, hi) = (a0.min(a1), a0.max(a1));
    let rate = |a: f64| -> f64 {
        let (c, d, dd) = rolling_centre(slant, radius, a);
        let x = Point::ORIGIN + slant.along * rolling_crossing(slant, radius, a);
        let touch = c - ogeom_math::Vector::Z * radius;
        let foot = c - u_v_reach(slant, radius, a);
        let triangle = (touch - x).cross(foot - x).magnitude() / 2.0;
        let tri_moment = ((x - Point::ORIGIN) + (touch - Point::ORIGIN) + (foot - Point::ORIGIN))
            * (triangle / 3.0);
        let angle = ((touch - c).dot(foot - c) / (radius * radius))
            .clamp(-1.0, 1.0)
            .acos();
        let segment = radius * radius * (angle - angle.sin()) / 2.0;
        let chord_mid = touch.midpoint(foot);
        let toward = (chord_mid - c) * (1.0 / chord_mid.distance(c));
        let lever = 4.0 * radius * (angle / 2.0).sin().powi(3) / (3.0 * (angle - angle.sin()));
        let seg_moment = ((c - Point::ORIGIN) + toward * lever) * segment;
        let area = triangle - segment;
        let moment = tri_moment - seg_moment;
        let speed = d.magnitude();
        let normal = d * (1.0 / speed);
        let turning = (dd - normal * dd.dot(normal)) * (1.0 / speed);
        speed * area - turning.dot(moment - (c - Point::ORIGIN) * area)
    };
    let nodes = gauss_legendre(40);
    let (m, h) = ((lo + hi) / 2.0, (hi - lo) / 2.0);
    nodes.iter().map(|&(x, w)| w * h * rate(m + h * x)).sum()
}

/// From the ball's centre to its contact on the block's top edge.
fn u_v_reach(slant: &Slant, radius: f64, a: f64) -> ogeom_math::Vector {
    let (s, c) = a.sin_cos();
    (slant.side * c + slant.top * s) * radius
}

/// How far `p` stands off the envelope of the rolling balls: its distance
/// to the nearest centre less the radius.
fn off_rolling_envelope(slant: &Slant, radius: f64, p: Point) -> f64 {
    let (a0, a1) = rolling_span(slant, radius);
    let (lo, hi) = (a0.min(a1) - 0.05, a0.max(a1) + 0.05);
    let dist = |a: f64| p.distance(rolling_centre(slant, radius, a).0);
    let best = (0..=400)
        .map(|i| lo + (hi - lo) * f64::from(i) / 400.0)
        .min_by(|x, y| dist(*x).total_cmp(&dist(*y)))
        .unwrap();
    let step = (hi - lo) / 400.0;
    let (mut a, mut b) = (best - step, best + step);
    let g = 0.5 * (5.0_f64.sqrt() - 1.0);
    for _ in 0..100 {
        let (x1, x2) = (b - g * (b - a), a + g * (b - a));
        if dist(x1) < dist(x2) {
            b = x2;
        } else {
            a = x1;
        }
    }
    dist((a + b) / 2.0) - radius
}

/// The fitted round among the faces a blend generated: how far its
/// surface stands off the envelope `off` measures, at most, over a grid
/// across it, and its area.
fn fitted_round(
    model: &ogeom_topo::Model,
    faces: &[ogeom_topo::Shape],
    off: &dyn Fn(Point) -> f64,
) -> (f64, f64) {
    use ogeom_geom::Surface as _;
    let mut worst: f64 = 0.0;
    let mut area = 0.0;
    for face in faces {
        let data = model.node(face).unwrap().data().as_face().unwrap();
        let surface = model.geometry().surface(data.surface).unwrap();
        if !matches!(surface, ogeom_geom::SurfaceGeometry::BSpline(_)) {
            continue;
        }
        let ((u0, u1), (v0, v1)) = surface.domain();
        for i in 0..=40 {
            for j in 0..=40 {
                let p = surface
                    .point_at(
                        u0 + (u1 - u0) * f64::from(i) / 40.0,
                        v0 + (v1 - v0) * f64::from(j) / 40.0,
                        T,
                    )
                    .unwrap();
                worst = worst.max(off(p).abs());
            }
        }
        area += ogeom_algo::surface_properties(model, face, fine(), T)
            .unwrap()
            .mass;
    }
    (worst, area)
}

/// A seam whose narrow strip is narrower than the ball's setback all
/// along: the ball touches the plate's top and rests on the strip's far
/// edge, the block's top edge, its centre on the ellipse where the plane a
/// radius over the top meets the cylinder a radius about the edge. The
/// fill is the volume between the seam and the balls, integrated over
/// their sections from the boxes' own planes and edge, and the round
/// stands a radius off the top and off the edge's line along its rails.
#[test]
fn a_round_rolls_on_the_far_edge_of_a_strip_narrower_than_its_setback() {
    let radius = 0.5;
    for alone in [true, false] {
        let mut model = ogeom_topo::Model::new();
        let slant = slanted_seam(&mut model, &SLANT);
        let before = volume(&model, &slant.shape);
        let rounded = if alone {
            ogeom_fillet::fillet_edge(&mut model, &slant.shape, &slant.seam, radius, T)
        } else {
            ogeom_fillet::fillet_edges(
                &mut model,
                &slant.shape,
                std::slice::from_ref(&slant.seam),
                radius,
                T,
            )
        }
        .unwrap();
        let diagnosis = ogeom_algo::check(&model, &rounded.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        // The round: the band fitted through the balls' arcs, within the
        // band's fit target of their envelope, its rails on the plate's top
        // and on the edge's line.
        let generated = rounded.history.generated(&slant.seam).to_vec();
        let (worst, area) = fitted_round(&model, &generated, &|p| {
            off_rolling_envelope(&slant, radius, p)
        });
        let (mut on_top, mut on_edge) = (false, false);
        for face in &generated {
            let data = model.node(face).unwrap().data().as_face().unwrap();
            if !matches!(
                model.geometry().surface(data.surface).unwrap(),
                ogeom_geom::SurfaceGeometry::BSpline(_)
            ) {
                continue;
            }
            for edge in explore(&model, face, Filter::OfType(ShapeType::Edge)).unwrap() {
                use ogeom_geom::Curve3d as _;
                let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = model
                    .node(&edge)
                    .unwrap()
                    .data()
                    .as_edge()
                    .unwrap()
                    .curve3d()
                else {
                    continue;
                };
                let curve = model.geometry().curve(*curve).unwrap();
                let points: Vec<Point> = (0..=32)
                    .map(|k| {
                        curve
                            .point_at(range.0 + (range.1 - range.0) * f64::from(k) / 32.0, T)
                            .unwrap()
                    })
                    .collect();
                // A rail set on the top's plane or the edge's line stands on
                // it to rounding: a nanometre is far below any fit error.
                let top = points.iter().all(|p| p.z.abs() < 1e-9);
                let line = points.iter().all(|p| {
                    let w = *p - Point::ORIGIN;
                    (w - slant.far * w.dot(slant.far)).magnitude() < 1e-9
                });
                // The band's fit target, 2e-4 at unit scale, is what a rail
                // is held to as well as the band.
                if top || line {
                    for p in &points {
                        assert!(
                            off_rolling_envelope(&slant, radius, *p).abs() < 2e-4,
                            "a rail point {p:?} stands off the balls"
                        );
                    }
                }
                on_top |= top;
                on_edge |= line;
            }
        }
        assert!(worst < 2e-4, "the round stands {worst} off the balls");
        assert!(
            on_top && on_edge,
            "the round's rails: top {on_top}, edge {on_edge}"
        );
        let want = rolling_fill(&slant, radius);
        let got = volume(&model, &rounded.shape) - before;
        // The band moves the fill by no more than its furthest stand-off
        // from the envelope over its area; the integral of the sections is
        // exact to rounding.
        assert!(
            (got - want).abs() <= worst * area,
            "the round added {got} against {want}, off {} with the band within {worst} over \
             {area}",
            got - want
        );
    }
}

/// Where the rolling run does not hold, the round is refused by name: the
/// seam running back to where the strip narrows to nothing, so the ball
/// comes onto the block's top beyond the strip's far edge and rounds the
/// corner between the top and the plate, and a post standing within the
/// ball's reach, beside a run resting on the edge all along and beside one
/// carried onto the strip.
#[test]
fn a_round_rolling_on_a_far_edge_refuses_where_the_ball_leaves_it() {
    for (seat, why) in [
        (
            Seat {
                ends: (0.0, 1.3),
                ..SLANT
            },
            "comes onto the face beyond that edge along the crease, rounding the corner between \
             that face and the one it touches rather than this edge",
        ),
        (
            Seat {
                post: true,
                ..SLANT
            },
            "comes within its radius of another face",
        ),
        (
            Seat {
                ends: (0.6, 1.9),
                turns: (-20.0, 20.0),
                post: true,
                ..SLANT
            },
            "comes within its radius of another face",
        ),
    ] {
        let ends = seat.ends;
        let mut model = ogeom_topo::Model::new();
        let slant = slanted_seam(&mut model, &seat);
        let refusal = ogeom_fillet::fillet_edge(&mut model, &slant.shape, &slant.seam, 0.5, T)
            .err()
            .unwrap_or_else(|| panic!("{ends:?}: the round was built"))
            .to_string();
        assert!(
            refusal.contains("rests on the far edge") && refusal.contains(why),
            "{ends:?}: {refusal}"
        );
    }
}

/// The strip's far edge parallel to the plate's top: the ball's centre
/// runs on a line, a radius over the top and a radius off the edge, and
/// the fill is a prism whose section, square to the seam, is the triangle
/// from the seam to the two contacts less the segment the ball's arc cuts
/// off it.
#[test]
fn a_round_rolls_on_a_far_edge_parallel_to_the_face_it_touches() {
    let (radius, lift) = (0.5, 0.3);
    let mut model = ogeom_topo::Model::new();
    let slant = slanted_seam(
        &mut model,
        &Seat {
            rise: 0.0,
            lift,
            turns: (0.0, 0.0),
            ..SLANT
        },
    );
    let before = volume(&model, &slant.shape);
    let rounded =
        ogeom_fillet::fillet_edge(&mut model, &slant.shape, &slant.seam, radius, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &rounded.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    // In the section x = 0: the edge at (0, lift), the seam where the side
    // face's plane through it meets z = 0, the centre a radius over z = 0
    // and a radius from the edge, on the side face's outward side.
    let (sy, sz) = (slant.side.y, slant.side.z);
    let seam = (sz * lift / sy, 0.0);
    let edge = (slant.lift.y, slant.lift.z);
    let centre = (-(radius * radius - (radius - lift).powi(2)).sqrt(), radius);
    let touch = (centre.0, 0.0);
    let triangle =
        ((touch.0 - seam.0) * (edge.1 - seam.1) - (touch.1 - seam.1) * (edge.0 - seam.0)).abs()
            / 2.0;
    let to = |p: (f64, f64)| (p.0 - centre.0, p.1 - centre.1);
    let (a, b) = (to(touch), to(edge));
    let angle = ((a.0 * b.0 + a.1 * b.1) / (radius * radius)).acos();
    let segment = radius * radius * (angle - angle.sin()) / 2.0;
    let want = (triangle - segment) * (slant.ends.1 - slant.ends.0);
    let generated = rounded.history.generated(&slant.seam).to_vec();
    let (worst, area) = fitted_round(&model, &generated, &|p| {
        (p.y - centre.0).hypot(p.z - centre.1) - radius
    });
    // The band's fit target, and the fill moved by no more than the band's
    // furthest stand-off over its area.
    assert!(worst < 2e-4, "the round stands {worst} off the balls");
    let got = volume(&model, &rounded.shape) - before;
    assert!(
        (got - want).abs() <= worst * area,
        "the round added {got} against {want}, with the band within {worst} over {area}"
    );
}

struct Tip {
    shape: ogeom_topo::Shape,
    seam: ogeom_topo::Shape,
    apex: Point,
    edges: [ogeom_math::Vector; 3],
}

/// A plate whose top is `z = 0`, fused with a box whose corner stands `h`
/// over it at the origin: the box's face spanned by its first two edges
/// out of the corner leans `lean` degrees off upright, those edges leaving
/// it at `turn` degrees and `turn + 90` round from level, and its third
/// edge runs into the box square to that face. The seam is where that face
/// meets the plate's top, the side of a right triangle whose right angle is
/// the corner.
fn box_tip(model: &mut ogeom_topo::Model, h: f64, lean: f64, turn: f64) -> Tip {
    use ogeom_math::{Direction, Vector};
    let (sa, ca) = lean.to_radians().sin_cos();
    let into = Vector::new(0.0, -ca, -sa);
    let level = Vector::X;
    let down = Vector::new(0.0, sa, -ca);
    let along = |deg: f64| {
        let (s, c) = deg.to_radians().sin_cos();
        level * c + down * s
    };
    let (d1, d2) = (along(turn), along(turn + 90.0));
    let apex = Point::new(0.0, 0.0, h);
    // A right-handed frame along the three edges, `into` its z.
    let x = if d1.cross(d2).dot(into) > 0.0 { d1 } else { d2 };
    let block = ogeom_algo::make_box(
        model,
        Frame::new(
            apex,
            Direction::new(into, T).unwrap(),
            Direction::new(x, T).unwrap(),
            T,
        )
        .unwrap(),
        (2.0, 2.0, 2.0),
        T,
    )
    .unwrap();
    let plate = ogeom_algo::make_box(
        model,
        Frame::new(Point::new(-3.0, -3.0, -3.0), Direction::Z, Direction::X, T).unwrap(),
        (6.0, 6.0, 3.0),
        T,
    )
    .unwrap();
    let shape = ogeom_bool::fuse(model, &block.shape, &plate.shape, T)
        .unwrap()
        .shape;
    let foot = |d: Vector| apex + d * (h / -d.z);
    let middle = foot(d1).midpoint(foot(d2));
    let seam = explore(model, &shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(model, e)
                .unwrap()
                .is_some_and(|(a, b)| {
                    let p = |v: &ogeom_topo::Shape| {
                        model
                            .node(v)
                            .and_then(|n| n.data().as_vertex().map(|d| d.point))
                            .unwrap()
                    };
                    p(&a).midpoint(p(&b)).distance(middle) < 1e-9
                })
        })
        .expect("the fused part has the seam");
    Tip {
        shape,
        seam,
        apex,
        edges: [d1, d2, into],
    }
}

/// One piece of a ball's run over the plate's top `z = 0`: its centre and
/// the centre's first and second derivatives at a parameter, the ball's
/// contact on what it holds off the plate for a centre, and the span of
/// the parameter it runs over.
struct Run<'a> {
    centre: &'a dyn Fn(f64) -> (Point, ogeom_math::Vector, ogeom_math::Vector),
    contact: &'a dyn Fn(Point) -> Point,
    span: (f64, f64),
}

/// The root of `f` in `[lo, hi]`, where it changes sign there.
fn bisect(f: &dyn Fn(f64) -> f64, lo: f64, hi: f64) -> f64 {
    let (mut a, mut b) = (lo, hi);
    let fa = f(a);
    for _ in 0..200 {
        let m = (a + b) / 2.0;
        if (f(m) > 0.0) == (fa > 0.0) {
            a = m;
        } else {
            b = m;
        }
    }
    (a + b) / 2.0
}

/// Where along a seam through `origin` along `along` the plane square to
/// a run's path at `s` crosses it.
fn crossing_of(run: &Run<'_>, origin: Point, along: ogeom_math::Vector, s: f64) -> f64 {
    let (c, d, _) = (run.centre)(s);
    d.dot(c - origin) / d.dot(along)
}

/// The volume between a seam on the plate's top, through `origin` along
/// `along`, and the balls of `radius` over `runs`, swept by the balls'
/// sections: each, in the plane square to the centre's path, is the
/// triangle from the seam to the two contacts less the circular segment
/// the ball's arc cuts off it, and a plane turning with the path sweeps
/// `(|c'| - n'.(p - c)) dA ds` at each point `p` of it.
fn swept_runs(runs: &[Run<'_>], origin: Point, along: ogeom_math::Vector, radius: f64) -> f64 {
    let nodes = gauss_legendre(40);
    let mut total = 0.0;
    for run in runs {
        let rate = |s: f64| -> f64 {
            let (c, d, dd) = (run.centre)(s);
            let x = origin + along * crossing_of(run, origin, along, s);
            let touch = c - ogeom_math::Vector::Z * radius;
            let foot = (run.contact)(c);
            let triangle = (touch - x).cross(foot - x).magnitude() / 2.0;
            let tri_moment =
                ((x - Point::ORIGIN) + (touch - Point::ORIGIN) + (foot - Point::ORIGIN))
                    * (triangle / 3.0);
            let angle = ((touch - c).dot(foot - c) / (radius * radius))
                .clamp(-1.0, 1.0)
                .acos();
            let segment = radius * radius * (angle - angle.sin()) / 2.0;
            let chord_mid = touch.midpoint(foot);
            let toward = (chord_mid - c) * (1.0 / chord_mid.distance(c));
            let lever = 4.0 * radius * (angle / 2.0).sin().powi(3) / (3.0 * (angle - angle.sin()));
            let seg_moment = ((c - Point::ORIGIN) + toward * lever) * segment;
            let area = triangle - segment;
            let moment = tri_moment - seg_moment;
            let speed = d.magnitude();
            let normal = d * (1.0 / speed);
            let turning = (dd - normal * dd.dot(normal)) * (1.0 / speed);
            speed * area - turning.dot(moment - (c - Point::ORIGIN) * area)
        };
        // Each piece's span runs either way of its parameter: its volume
        // is the integral's size.
        let (lo, hi) = run.span;
        let pieces = 16;
        let mut piece = 0.0;
        for k in 0..pieces {
            let a = lo + (hi - lo) * f64::from(k) / f64::from(pieces);
            let b = lo + (hi - lo) * f64::from(k + 1) / f64::from(pieces);
            let (m, h) = ((a + b) / 2.0, (b - a) / 2.0);
            piece += nodes
                .iter()
                .map(|&(x, w)| w * h * rate(m + h * x))
                .sum::<f64>();
        }
        total += piece.abs();
    }
    total
}

/// How far `p` stands off the envelope of the balls of `radius` over
/// `runs`: its distance to the nearest centre less the radius.
fn off_runs(runs: &[Run<'_>], radius: f64, p: Point) -> f64 {
    runs.iter()
        .map(|run| {
            let (lo, hi) = run.span;
            let dist = |s: f64| p.distance((run.centre)(s).0);
            let at = |i: i32| lo + (hi - lo) * f64::from(i) / 400.0;
            let best = (0..=400)
                .min_by(|x, y| dist(at(*x)).total_cmp(&dist(at(*y))))
                .unwrap();
            let (mut a, mut b) = (at((best - 1).max(0)), at((best + 1).min(400)));
            let g = 0.5 * (5.0_f64.sqrt() - 1.0);
            for _ in 0..100 {
                let (x1, x2) = (b - g * (b - a), a + g * (b - a));
                if dist(x1) < dist(x2) {
                    b = x2;
                } else {
                    a = x1;
                }
            }
            dist((a + b) / 2.0)
        })
        .fold(f64::INFINITY, f64::min)
        - radius
}

/// Round `seam` alone: the result valid, its fitted bands within the
/// band's fit target of the envelope of `runs`' balls, and the volume it
/// added equal to their swept sections within what the bands' stand-off
/// over their area can move it. Returns how many fitted bands it has.
fn holds_to_runs(
    model: &mut ogeom_topo::Model,
    shape: &ogeom_topo::Shape,
    seam: &ogeom_topo::Shape,
    radius: f64,
    runs: &[Run<'_>],
    (origin, along): (Point, ogeom_math::Vector),
) -> usize {
    let before = volume(model, shape);
    let rounded = ogeom_fillet::fillet_edge(model, shape, seam, radius, T).unwrap();
    let diagnosis = ogeom_algo::check(model, &rounded.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let generated = rounded.history.generated(seam).to_vec();
    let bands = generated
        .iter()
        .filter(|face| {
            let data = model.node(face).unwrap().data().as_face().unwrap();
            matches!(
                model.geometry().surface(data.surface).unwrap(),
                ogeom_geom::SurfaceGeometry::BSpline(_)
            )
        })
        .count();
    let (worst, area) = fitted_round(model, &generated, &|p| off_runs(runs, radius, p));
    // The band's fit target, 2e-4 at unit scale.
    assert!(worst < 2e-4, "the round stands {worst} off the balls");
    let want = swept_runs(runs, origin, along, radius);
    let got = volume(model, &rounded.shape) - before;
    assert!(
        (got - want).abs() <= worst * area,
        "the round added {got} against {want}, off {} with the bands within {worst} over \
         {area}",
        got - want
    );
    bands
}

/// The seam run on to where the strip beside it is wider than the ball's
/// setback: the ball rests on the strip's far edge, its centre on the
/// ellipse, until it comes down on the strip where the ellipse meets the
/// line a radius over the plate and a radius off the strip, tangent to it,
/// and runs on that line, the cylinder between the two planes, to the
/// seam's end. The fill is the two pieces' swept sections, from the boxes'
/// own planes and edge. The round is two bands meeting on the section
/// where the ball holds all three.
#[test]
fn a_round_rolls_onto_a_strip_that_widens_past_its_setback() {
    use ogeom_math::Vector;
    let radius = 0.5;
    let mut model = ogeom_topo::Model::new();
    let slant = slanted_seam(
        &mut model,
        &Seat {
            ends: (0.6, 1.9),
            turns: (-20.0, 20.0),
            ..SLANT
        },
    );
    let rolled = |a: f64| rolling_centre(&slant, radius, a);
    let on_edge = |c: Point| slant.lift + slant.far * (c - slant.lift).dot(slant.far);
    // The ball comes down on the strip where its angle round the edge
    // reaches the strip's normal.
    let switch = rolling_crossing(&slant, radius, 0.0);
    assert!(
        switch > slant.ends.0 && switch < slant.ends.1,
        "the ball comes down at {switch}"
    );
    let start = bisect(
        &|a| rolling_crossing(&slant, radius, a) - slant.ends.0,
        0.0,
        core::f64::consts::FRAC_PI_2,
    );
    // A radius over the plate and a radius off the strip, square to the
    // seam.
    let n = slant.side;
    let w = {
        let m = Vector::Z.dot(n);
        let (a, b) = (
            (radius - m * radius) / (1.0 - m * m),
            (radius - m * radius) / (1.0 - m * m),
        );
        Vector::Z * a + n * b
    };
    let on_line = |x: f64| {
        (
            Point::ORIGIN + slant.along * x + w,
            slant.along,
            Vector::ZERO,
        )
    };
    assert!(
        rolled(0.0).0.distance(on_line(switch).0) < 1e-9,
        "the ellipse meets the line where the ball comes down"
    );
    let on_strip = |c: Point| c - n * radius;
    let runs = [
        Run {
            centre: &rolled,
            contact: &on_edge,
            span: (start, 0.0),
        },
        Run {
            centre: &on_line,
            contact: &on_strip,
            span: (switch, slant.ends.1),
        },
    ];
    let bands = holds_to_runs(
        &mut model,
        &slant.shape,
        &slant.seam,
        radius,
        &runs,
        (Point::ORIGIN, slant.along),
    );
    assert_eq!(bands, 2, "a band for each piece");
}

/// A box's corner standing over a plate, its steep face a right triangle
/// over the seam narrower than the ball's setback. From the seam's far
/// end the ball touches the plate and rests on the triangle's second edge
/// out of the corner, its centre on the ellipse, until its contact on the
/// edge reaches the corner. It then rolls round the corner to the seam's
/// near end, its centre on the circle where the plane a radius over the
/// plate meets the sphere a radius about the corner, tangent there to the
/// ellipse. The fill is the two pieces' swept sections, from the box's own
/// corner, edge and faces, and the corner's piece is a patch of a spindle
/// torus whose arcs all end on the corner.
#[test]
fn a_round_rolls_round_the_corner_of_a_narrow_triangle() {
    use ogeom_math::Vector;
    let (radius, h) = (0.38, 0.1);
    let mut model = ogeom_topo::Model::new();
    let tip = box_tip(&mut model, h, 4.0, 51.0);
    let [d1, d2, into] = tip.edges;
    let w = tip.apex;
    let foot = |d: Vector| w + d * (h / -d.z);
    let (p1, p2) = (foot(d1), foot(d2));
    let along = (p2 - p1) * (1.0 / p1.distance(p2));
    let length = p1.distance(p2);
    // The edge out of the corner along `d2`, between the triangle and the
    // box's face square to `d1`: the ball a radius over the plate resting
    // on it at angle `a` round it from the triangle's outward normal `u`
    // toward that face's, `-d1`.
    let u = -into;
    let v = {
        let v = -d1 - u * (-d1).dot(u);
        v * (1.0 / v.magnitude())
    };
    let on_edge = |a: f64| {
        let (s, c) = a.sin_cos();
        let l = (radius - h - radius * (c * u.z + s * v.z)) / d2.z;
        let dl = radius * (s * u.z - c * v.z) / d2.z;
        let ddl = radius * (c * u.z + s * v.z) / d2.z;
        (
            w + d2 * l + (u * c + v * s) * radius,
            d2 * dl + (v * c - u * s) * radius,
            d2 * ddl - (u * c + v * s) * radius,
        )
    };
    let edge_foot = |c: Point| w + d2 * (c - w).dot(d2);
    // A radius over the plate and a radius from the corner.
    let rho = (h * (2.0 * radius - h)).sqrt();
    let o = Point::new(w.x, w.y, radius);
    let round = |t: f64| {
        let (s, c) = t.sin_cos();
        (
            o + Vector::new(c, s, 0.0) * rho,
            Vector::new(-s, c, 0.0) * rho,
            Vector::new(-c, -s, 0.0) * rho,
        )
    };
    let corner = |_: Point| w;
    let x_of = |f: &dyn Fn(f64) -> (Point, Vector, Vector), s: f64| {
        let (c, d, _) = f(s);
        d.dot(c - p1) / d.dot(along)
    };
    // The arc of angles on which the edge alone holds the ball, from the
    // triangle's normal to the next face's; the foot reaches the corner
    // where `l(a) = 0`, `cos a u.z + sin a v.z = (r - h) / r`.
    let phi = (-d1).dot(v).atan2((-d1).dot(u));
    let samples = |lo: f64, hi: f64| {
        (0..400).map(move |i| {
            (
                lo + (hi - lo) * f64::from(i) / 400.0,
                lo + (hi - lo) * f64::from(i + 1) / 400.0,
            )
        })
    };
    let lift = |a: f64| (on_edge(a).0 - w).dot(d2);
    let a_corner = samples(0.0, phi)
        .find(|(a, b)| (lift(*a) > 0.0) != (lift(*b) > 0.0))
        .map(|(a, b)| bisect(&lift, a, b))
        .expect("the foot reaches the corner");
    let to_end = |a: f64| x_of(&on_edge, a) - length;
    let a_end = samples(0.0, phi)
        .find(|(a, b)| (to_end(*a) > 0.0) != (to_end(*b) > 0.0))
        .map(|(a, b)| bisect(&to_end, a, b))
        .expect("the edge holds the ball at the seam's far end");
    let t_corner = {
        let c = on_edge(a_corner).0 - o;
        c.y.atan2(c.x)
    };
    assert!(
        round(t_corner).0.distance(on_edge(a_corner).0) < 1e-9,
        "the circle meets the ellipse where the foot reaches the corner"
    );
    // Round the corner the way the sections move toward the seam's near
    // end, to where they reach it.
    let toward = -(x_of(&round, t_corner + 1e-6) - x_of(&round, t_corner - 1e-6)).signum();
    let to_start = |t: f64| x_of(&round, t);
    let t_start = samples(t_corner, t_corner + toward * core::f64::consts::PI)
        .find(|(a, b)| (to_start(*a) > 0.0) != (to_start(*b) > 0.0))
        .map(|(a, b)| bisect(&to_start, a, b))
        .expect("the corner holds the ball at the seam's near end");
    let runs = [
        Run {
            centre: &on_edge,
            contact: &edge_foot,
            span: (a_end, a_corner),
        },
        Run {
            centre: &round,
            contact: &corner,
            span: (t_corner, t_start),
        },
    ];
    let bands = holds_to_runs(
        &mut model,
        &tip.shape,
        &tip.seam,
        radius,
        &runs,
        (p1, along),
    );
    assert_eq!(bands, 2, "a band for each piece");
}
