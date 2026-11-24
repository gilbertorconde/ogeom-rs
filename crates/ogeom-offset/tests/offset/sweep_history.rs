//! Each sweep's history, against the sweep's own definition: a face a
//! profile edge is recorded as generating is one whose every sampled point
//! lies on the path that edge sweeps out, and every such face is recorded.
//! The paths are written here from each operation's documented motion (a
//! screw, a translation turned at a mitre, a revolution), not read from
//! the construction.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_algo::History;
use ogeom_algo::sweep::roles;
use ogeom_core::{Role, Tolerances};
use ogeom_geom::{PlaneSurface, SurfaceGeometry};
use ogeom_math::{Axis, Direction, Plane, Point, Transform, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A motion of the profile: where its point `p` stands at sweep parameter
/// `s`. Affine in `p` for each `s`, so a straight profile edge stays a
/// segment.
type Motion = Box<dyn Fn(f64, Point) -> Point>;

/// One run of the sweep: a motion over a parameter range.
struct Run {
    motion: Motion,
    range: (f64, f64),
}

/// The distance from `p` to segment `a`-`b`.
fn to_segment(p: Point, a: Point, b: Point) -> f64 {
    let d = b - a;
    let t = ((p - a).dot(d) / d.dot(d)).clamp(0.0, 1.0);
    p.distance(a + d * t)
}

/// The distance from `p` to the surface the segment `a`-`b` sweeps under
/// `runs`: a dense scan of each run's parameter, then a golden-section
/// search about every scanned dip.
fn to_swept(p: Point, a: Point, b: Point, runs: &[Run]) -> f64 {
    let mut best = f64::INFINITY;
    for run in runs {
        let at = |s: f64| to_segment(p, (run.motion)(s, a), (run.motion)(s, b));
        let n = 240;
        let (lo, hi) = run.range;
        let step = (hi - lo) / f64::from(n);
        let scan: Vec<f64> = (0..=n).map(|k| at(lo + step * f64::from(k))).collect();
        for k in 0..=n as usize {
            let dips =
                (k == 0 || scan[k] <= scan[k - 1]) && (k == n as usize || scan[k] <= scan[k + 1]);
            if !dips {
                continue;
            }
            #[allow(clippy::cast_precision_loss)]
            let centre = lo + step * k as f64;
            let (mut x0, mut x1) = ((centre - step).max(lo), (centre + step).min(hi));
            let g = 0.5 * (5.0_f64.sqrt() - 1.0);
            for _ in 0..80 {
                let (c, d) = (x1 - g * (x1 - x0), x0 + g * (x1 - x0));
                if at(c) < at(d) {
                    x1 = d;
                } else {
                    x0 = c;
                }
            }
            best = best.min(at(0.5 * (x0 + x1))).min(scan[k]);
        }
    }
    best
}

/// A plane face through the corners, in order.
fn polygon_face(model: &mut Model, corners: &[Point]) -> Shape {
    let wire = ogeom_algo::make_polygon(model, corners, true, T)
        .unwrap()
        .shape;
    let plane = ogeom_algo::find_plane(model, &wire, T).unwrap().unwrap();
    let surface: SurfaceGeometry = PlaneSurface::new(plane).into();
    ogeom_algo::make_face(model, surface, &[wire], T)
        .unwrap()
        .shape
}

/// A straight wire through the points.
fn path(model: &mut Model, points: &[Point]) -> Shape {
    ogeom_algo::make_polygon(model, points, false, T)
        .unwrap()
        .shape
}

/// The two ends of a straight profile edge, placed.
fn ends_of(model: &Model, edge: &Shape) -> (Point, Point) {
    let (curve, range) = ogeom_algo::placed_edge_curve(model, edge, T)
        .unwrap()
        .unwrap();
    use ogeom_geom::Curve3d as _;
    (
        curve.point_at(range.0, T).unwrap(),
        curve.point_at(range.1, T).unwrap(),
    )
}

/// Points of `face`: its mesh's nodes, and each triangle's centre carried
/// onto the face's surface, so a face meshed by its corners alone still
/// shows its inside.
fn samples(model: &Model, face: &Shape) -> Vec<Point> {
    let deflection = ogeom_mesh::Deflection {
        chord: 1e-2,
        ..ogeom_mesh::Deflection::default()
    };
    let mesh = ogeom_mesh::triangulate(model, face, deflection, T).unwrap();
    let data = model.node(face).unwrap().data().as_face().unwrap();
    let surface = model.geometry().surface(data.surface).unwrap();
    let placement = face.transform(model.datums()).unwrap();
    let inverse = placement.inverse().unwrap();
    let mut out = mesh.positions.clone();
    for tri in &mesh.triangles {
        let [a, b, c] = tri.map(|i| mesh.positions[i as usize]);
        let centre = Point::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        );
        let foot = ogeom_algo::project_on_surface(surface, inverse.apply(centre), 16, T).unwrap();
        out.push(placement.apply(foot.point));
    }
    // At most a few dozen, spread over the mesh's order.
    let step = out.len().div_ceil(48).max(1);
    out.into_iter().step_by(step).collect()
}

/// The role a shape's provenance gives it, and the entities it names.
fn role_of(model: &Model, shape: &Shape) -> Option<(Role, Vec<ogeom_core::EntityId>)> {
    model
        .provenance_of(shape)
        .and_then(|p| p.role().map(|r| (r, p.inputs().to_vec())))
}

/// What a sweep's history was measured to do.
struct Measured {
    /// The farthest a sampled point of a recorded wall stands off its
    /// edge's swept surface.
    worst_on: f64,
    /// The nearest any face not recorded for an edge comes to lying
    /// wholly on that edge's swept surface: its farthest sample from it.
    least_off: f64,
}

/// Checks the history of a sweep whose result is `shape`:
///
/// - each profile edge generates exactly the faces of `shape` whose every
///   sample lies within `bound` of its swept surface (`runs`), and at
///   least one;
/// - every face of `shape` is traced from the profile, one of its edges,
///   or `others`;
/// - the faces carrying [`roles::SWEEP_BOTTOM`] and [`roles::SWEEP_TOP`]
///   name `bottom_from` and `top_from` and lie in `bottom_plane` and
///   `top_plane`.
#[allow(clippy::too_many_arguments)]
fn check_sweep(
    model: &Model,
    shape: &Shape,
    history: &History,
    profile: &Shape,
    others: Option<&[&Shape]>,
    runs: &[Run],
    bound: f64,
    caps: &[(Role, &Shape, Option<Plane>)],
) -> Measured {
    let edges = explore_unique(model, profile, ShapeType::Edge).unwrap();
    let faces = explore_unique(model, shape, ShapeType::Face).unwrap();
    let mut measured = Measured {
        worst_on: 0.0,
        least_off: f64::INFINITY,
    };
    let face_samples: Vec<Vec<Point>> = faces.iter().map(|f| samples(model, f)).collect();
    for edge in &edges {
        let (a, b) = ends_of(model, edge);
        let generated = history.generated(edge);
        let mut walls = 0;
        for (face, points) in faces.iter().zip(&face_samples) {
            let off = points
                .iter()
                .map(|p| to_swept(*p, a, b, runs))
                .fold(0.0, f64::max);
            let recorded = generated.iter().any(|g| g.is_same(face));
            if recorded {
                walls += 1;
                measured.worst_on = measured.worst_on.max(off);
                if off > bound {
                    let others: Vec<f64> = edges
                        .iter()
                        .map(|o| {
                            let (a, b) = ends_of(model, o);
                            points
                                .iter()
                                .map(|p| to_swept(*p, a, b, runs))
                                .fold(0.0, f64::max)
                        })
                        .collect();
                    let median = {
                        let mut d: Vec<f64> =
                            points.iter().map(|p| to_swept(*p, a, b, runs)).collect();
                        d.sort_by(f64::total_cmp);
                        d[d.len() / 2]
                    };
                    panic!(
                        "a recorded wall stands {off} off its edge's swept surface (median {median}); against each edge {others:?}"
                    );
                }
            } else {
                measured.least_off = measured.least_off.min(off);
                assert!(
                    off > bound,
                    "a face lying on an edge's swept surface (within {off}) is not recorded: {:?}, role {:?}",
                    ogeom_algo::shape_bounds(model, face, T).unwrap(),
                    role_of(model, face)
                );
            }
        }
        assert!(walls > 0, "an edge generates no face");
    }
    // Every face is reached from some input.
    let mut inputs: Vec<&Shape> = vec![profile];
    inputs.extend(edges.iter());
    let vertices = explore_unique(model, profile, ShapeType::Vertex).unwrap();
    inputs.extend(vertices.iter());
    inputs.extend(others.unwrap_or_default().iter().copied());
    for face in faces.iter().filter(|_| others.is_some()) {
        let reached = inputs.iter().any(|i| {
            history.generated(i).iter().any(|g| g.is_same(face))
                || history.modified(i).iter().any(|g| g.is_same(face))
        });
        assert!(
            reached,
            "a face of the sweep is traced from no input: {:?}, role {:?}",
            ogeom_algo::shape_bounds(model, face, T).unwrap(),
            role_of(model, face)
        );
    }
    // The caps.
    for (role, from, plane) in caps {
        let with_role: Vec<(&Shape, &Vec<Point>)> = faces
            .iter()
            .zip(&face_samples)
            .filter(|(f, _)| role_of(model, f).is_some_and(|(r, _)| r == *role))
            .collect();
        assert!(!with_role.is_empty(), "no face carries the role {role:?}");
        for (face, points) in with_role {
            let (_, named) = role_of(model, face).unwrap();
            assert_eq!(
                named,
                model.identity_of(from).into_iter().collect::<Vec<_>>()
            );
            if let Some(plane) = plane {
                let off = points
                    .iter()
                    .map(|p| plane.distance_to(*p))
                    .fold(0.0, f64::max);
                assert!(off <= bound, "a cap stands {off} off its plane");
            }
        }
    }
    measured
}

/// A square of side `2 * half` about `centre`, in the plane of `u` and
/// `v`.
fn square_corners(centre: Point, u: Vector, v: Vector, half: f64) -> Vec<Point> {
    [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
        .iter()
        .map(|(a, b)| centre + u * (a * half) + v * (b * half))
        .collect()
}

#[test]
fn a_helical_sweep_traces_each_wall_to_its_edge() {
    let mut model = Model::new();
    // A unit square in the XZ plane, clear of the z axis, screwed one and
    // a half turns at a pitch of three.
    let corners = [
        Point::new(2.0, 0.0, 0.0),
        Point::new(3.0, 0.0, 0.0),
        Point::new(3.0, 0.0, 1.0),
        Point::new(2.0, 0.0, 1.0),
    ];
    let profile = polygon_face(&mut model, &corners);
    let axis = Axis::new(Point::ORIGIN, Direction::Z);
    let (pitch, turns) = (3.0, 1.5);
    let built =
        ogeom_offset::make_helical_sweep(&mut model, &profile, axis, pitch, turns, false, 0.0, T)
            .unwrap();
    let total = core::f64::consts::TAU * turns;
    let screw = move |s: f64, p: Point| -> Point {
        Transform::translation(Vector::new(0.0, 0.0, pitch * s / core::f64::consts::TAU))
            .apply(Transform::rotation(axis, s).apply(p))
    };
    let runs = [Run {
        motion: Box::new(screw),
        range: (0.0, total),
    }];
    let start = Plane::through(Point::ORIGIN, Direction::Y);
    let end = {
        let n = Vector::new(-total.sin(), total.cos(), 0.0);
        Plane::through(screw(total, Point::ORIGIN), Direction::new(n, T).unwrap())
    };
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &profile,
        Some(&[]),
        &runs,
        1e-3,
        &[
            (roles::SWEEP_BOTTOM, &profile, Some(start)),
            (roles::SWEEP_TOP, &profile, Some(end)),
        ],
    );
    eprintln!(
        "helical: recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

/// A run carrying the profile by `turn` and then along `dir` from `from`
/// (the point the profile's own `centre` is carried to), over `range`.
fn leg(turn: Transform, centre: Point, from: Point, dir: Vector, range: (f64, f64)) -> Run {
    Run {
        motion: Box::new(move |s, p| from + dir * s + (turn.apply(p) - turn.apply(centre))),
        range,
    }
}

/// A run turning the profile, already carried by `before`, about `axis`
/// through up to `angle`.
fn turn_run(before: Transform, axis: Axis, angle: f64) -> Run {
    Run {
        motion: Box::new(move |s, p| Transform::rotation(axis, s).apply(before.apply(p))),
        range: (0.0, angle),
    }
}

/// A square of half-width a half in the XY plane about the origin: the
/// section of the pipes below, square to a spine starting up the z axis.
fn xy_square(model: &mut Model) -> Shape {
    let corners = square_corners(
        Point::ORIGIN,
        Vector::new(1.0, 0.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
        0.5,
    );
    polygon_face(model, &corners)
}

#[test]
fn a_mitred_pipe_traces_each_piece_of_a_wall_to_its_edge() {
    // Up the z axis, then along x: under the rotation-minimizing frame the
    // section turns a quarter turn about y at the corner, and each wall is
    // its leg's prism run on to the plane halving the corner.
    let mut model = Model::new();
    let profile = xy_square(&mut model);
    let (up, along) = (3.0, 2.0);
    let corner = Point::new(0.0, 0.0, up);
    let spine = path(
        &mut model,
        &[Point::ORIGIN, corner, Point::new(along, 0.0, up)],
    );
    let built = ogeom_offset::make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &ogeom_offset::PipeLaw::RotationMinimizing,
        ogeom_offset::PipeCorners::Mitre,
        1e-4,
        T,
    )
    .unwrap();
    let quarter = Transform::rotation(
        Axis::new(Point::ORIGIN, Direction::Y),
        core::f64::consts::FRAC_PI_2,
    );
    let runs = [
        leg(
            Transform::IDENTITY,
            Point::ORIGIN,
            Point::ORIGIN,
            Vector::new(0.0, 0.0, 1.0),
            (0.0, up + 1.0),
        ),
        leg(
            quarter,
            Point::ORIGIN,
            corner,
            Vector::new(1.0, 0.0, 0.0),
            (-1.0, along),
        ),
    ];
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &profile,
        Some(&[&spine]),
        &runs,
        1e-6,
        &[
            (
                roles::SWEEP_BOTTOM,
                &profile,
                Some(Plane::through(Point::ORIGIN, Direction::Z)),
            ),
            (
                roles::SWEEP_TOP,
                &profile,
                Some(Plane::through(Point::new(along, 0.0, up), Direction::X)),
            ),
        ],
    );
    // Each side wall is two pieces, one down either leg; the walls in the
    // plane of the turn are one face each, their legs' prisms meeting in
    // one plane.
    let edges = explore_unique(&model, &profile, ShapeType::Edge).unwrap();
    let mut counts: Vec<usize> = edges
        .iter()
        .map(|e| {
            built
                .history
                .generated(e)
                .iter()
                .filter(|g| model.kind_of(g).unwrap() == ShapeType::Face)
                .count()
        })
        .collect();
    counts.sort_unstable();
    eprintln!(
        "mitred pipe: walls per edge {counts:?}, recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

#[test]
fn a_round_cornered_pipe_traces_its_corner_to_the_edges() {
    // As above, with the corner rounded: each leg ends square at the corner
    // and the section's half outside the turn revolves about the corner's
    // axis, a line along y through the corner offset outward... the axis
    // stands at the corner itself for a section turned about its centre.
    let mut model = Model::new();
    let profile = xy_square(&mut model);
    let (up, along) = (3.0, 2.0);
    let corner = Point::new(0.0, 0.0, up);
    let spine = path(
        &mut model,
        &[Point::ORIGIN, corner, Point::new(along, 0.0, up)],
    );
    let built = ogeom_offset::make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &ogeom_offset::PipeLaw::RotationMinimizing,
        ogeom_offset::PipeCorners::Round,
        1e-4,
        T,
    )
    .unwrap();
    let quarter = Transform::rotation(
        Axis::new(Point::ORIGIN, Direction::Y),
        core::f64::consts::FRAC_PI_2,
    );
    let lift = Transform::translation(corner.to_vector());
    let runs = [
        leg(
            Transform::IDENTITY,
            Point::ORIGIN,
            Point::ORIGIN,
            Vector::new(0.0, 0.0, 1.0),
            (0.0, up),
        ),
        turn_run(
            lift,
            Axis::new(corner, Direction::Y),
            core::f64::consts::FRAC_PI_2,
        ),
        leg(
            quarter,
            Point::ORIGIN,
            corner,
            Vector::new(1.0, 0.0, 0.0),
            (0.0, along),
        ),
    ];
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &profile,
        Some(&[&spine]),
        &runs,
        1e-6,
        &[
            (
                roles::SWEEP_BOTTOM,
                &profile,
                Some(Plane::through(Point::ORIGIN, Direction::Z)),
            ),
            (
                roles::SWEEP_TOP,
                &profile,
                Some(Plane::through(Point::new(along, 0.0, up), Direction::X)),
            ),
        ],
    );
    eprintln!(
        "round-cornered pipe: recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

#[test]
fn a_pipe_round_a_closed_spine_traces_every_wall() {
    // A square loop in the XY plane started halfway along its first side,
    // the section square to it: no caps, each wall mitred at four corners,
    // the section a quarter turn further round about z down each leg.
    let mut model = Model::new();
    let c = Point::new(2.0, 0.0, 0.0);
    let corners = square_corners(
        c,
        Vector::new(0.0, 1.0, 0.0),
        Vector::new(0.0, 0.0, 1.0),
        0.5,
    );
    let profile = polygon_face(&mut model, &corners);
    let loop_points = [
        c,
        Point::new(4.0, 0.0, 0.0),
        Point::new(4.0, 4.0, 0.0),
        Point::new(0.0, 4.0, 0.0),
        Point::ORIGIN,
        c,
    ];
    let spine = ogeom_algo::make_polygon(&mut model, &loop_points[..5], true, T)
        .unwrap()
        .shape;
    let built = ogeom_offset::make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &ogeom_offset::PipeLaw::RotationMinimizing,
        ogeom_offset::PipeCorners::Mitre,
        1e-4,
        T,
    )
    .unwrap();
    let dirs = [
        Vector::new(1.0, 0.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
        Vector::new(-1.0, 0.0, 0.0),
        Vector::new(0.0, -1.0, 0.0),
        Vector::new(1.0, 0.0, 0.0),
    ];
    let runs: Vec<Run> = (0..5)
        .map(|k| {
            #[allow(clippy::cast_precision_loss)]
            let turn = Transform::rotation(
                Axis::new(c, Direction::Z),
                core::f64::consts::FRAC_PI_2 * k as f64,
            );
            let length = loop_points[k].distance(loop_points[k + 1]);
            leg(turn, c, loop_points[k], dirs[k], (-1.0, length + 1.0))
        })
        .collect();
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &profile,
        Some(&[&spine]),
        &runs,
        1e-6,
        &[],
    );
    eprintln!(
        "closed pipe: recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

#[test]
fn a_pipe_under_scale_and_twist_laws_traces_each_wall() {
    // Straight up the z axis, the section shrinking to half and turning a
    // quarter turn counter-clockwise about the travel, both linearly in the
    // length run.
    let mut model = Model::new();
    let profile = xy_square(&mut model);
    let length = 4.0;
    let spine = path(&mut model, &[Point::ORIGIN, Point::new(0.0, 0.0, length)]);
    let laws = ogeom_offset::PipeLaws {
        scale: ogeom_offset::SpineLaw::Linear(0.5),
        twist: ogeom_offset::SpineLaw::Linear(core::f64::consts::FRAC_PI_2),
        ..ogeom_offset::PipeLaws::default()
    };
    let tolerance = 1e-4;
    let built = ogeom_offset::make_pipe_shell_laws(
        &mut model,
        &profile,
        &spine,
        &ogeom_offset::PipeLaw::RotationMinimizing,
        &laws,
        ogeom_offset::PipeCorners::Mitre,
        tolerance,
        T,
    )
    .unwrap();
    let runs = [Run {
        motion: Box::new(move |s, p| {
            let share = s / length;
            let scaled = Point::new(p.x * (1.0 - 0.5 * share), p.y * (1.0 - 0.5 * share), 0.0);
            let turned = Transform::rotation(
                Axis::new(Point::ORIGIN, Direction::Z),
                core::f64::consts::FRAC_PI_2 * share,
            )
            .apply(scaled);
            turned + Vector::new(0.0, 0.0, s)
        }),
        range: (0.0, length),
    }];
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &profile,
        Some(&[&spine]),
        &runs,
        tolerance * 2.0,
        &[
            (
                roles::SWEEP_BOTTOM,
                &profile,
                Some(Plane::through(Point::ORIGIN, Direction::Z)),
            ),
            (
                roles::SWEEP_TOP,
                &profile,
                Some(Plane::through(Point::new(0.0, 0.0, length), Direction::Z)),
            ),
        ],
    );
    eprintln!(
        "law pipe: recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

#[test]
fn a_revolution_up_to_a_face_traces_its_walls_through_the_common() {
    // The unit square of the XZ plane clear of the z axis, turned about it
    // until it meets the plane through the axis at sixty degrees: that
    // plane's face is the far cap.
    let mut model = Model::new();
    let corners = [
        Point::new(2.0, 0.0, 0.0),
        Point::new(3.0, 0.0, 0.0),
        Point::new(3.0, 0.0, 1.0),
        Point::new(2.0, 0.0, 1.0),
    ];
    let profile = polygon_face(&mut model, &corners);
    let stop = core::f64::consts::FRAC_PI_3;
    let radial = Vector::new(stop.cos(), stop.sin(), 0.0);
    let limit = polygon_face(
        &mut model,
        &[
            Point::ORIGIN + radial * 1.0 - Vector::new(0.0, 0.0, 2.0),
            Point::ORIGIN + radial * 5.0 - Vector::new(0.0, 0.0, 2.0),
            Point::ORIGIN + radial * 5.0 + Vector::new(0.0, 0.0, 3.0),
            Point::ORIGIN + radial * 1.0 + Vector::new(0.0, 0.0, 3.0),
        ],
    );
    let axis = Axis::new(Point::ORIGIN, Direction::Z);
    let built = ogeom_offset::make_revolution_until(&mut model, &profile, axis, &limit, T).unwrap();
    let runs = [Run {
        motion: Box::new(move |s, p| Transform::rotation(axis, s).apply(p)),
        range: (0.0, stop),
    }];
    let limit_plane = Plane::through(
        Point::ORIGIN,
        Direction::new(Vector::new(-stop.sin(), stop.cos(), 0.0), T).unwrap(),
    );
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &profile,
        Some(&[&limit]),
        &runs,
        1e-6,
        &[
            (
                roles::SWEEP_BOTTOM,
                &profile,
                Some(Plane::through(Point::ORIGIN, Direction::Y)),
            ),
            (roles::SWEEP_TOP, &limit, Some(limit_plane)),
        ],
    );
    eprintln!(
        "revolution up to a face: recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

#[test]
fn an_evolved_sweep_traces_each_wall_and_its_corner_to_the_edge() {
    // A moulding: a unit square two out from an L spine in the XY plane,
    // standing in the YZ plane outside the turn. Each wall is a prism down
    // the first leg, the square turned about the corner, and a prism down
    // the second leg.
    let mut model = Model::new();
    let (first, second) = (3.0, 2.0);
    let spine = path(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(first, 0.0, 0.0),
            Point::new(first, second, 0.0),
        ],
    );
    let profile = ogeom_algo::make_polygon(
        &mut model,
        &[
            Point::new(0.0, -2.0, 0.0),
            Point::new(0.0, -1.0, 0.0),
            Point::new(0.0, -1.0, 1.0),
            Point::new(0.0, -2.0, 1.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    let built = ogeom_offset::make_evolved(&mut model, &spine, &profile, T).unwrap();
    let corner = Point::new(first, 0.0, 0.0);
    let quarter = Transform::rotation(
        Axis::new(Point::ORIGIN, Direction::Z),
        core::f64::consts::FRAC_PI_2,
    );
    let runs = [
        leg(
            Transform::IDENTITY,
            Point::ORIGIN,
            Point::ORIGIN,
            Vector::new(1.0, 0.0, 0.0),
            (0.0, first),
        ),
        turn_run(
            Transform::translation(corner.to_vector()),
            Axis::new(corner, Direction::Z),
            core::f64::consts::FRAC_PI_2,
        ),
        leg(
            quarter,
            Point::ORIGIN,
            corner,
            Vector::new(0.0, 1.0, 0.0),
            (0.0, second),
        ),
    ];
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &profile,
        Some(&[&spine]),
        &runs,
        1e-6,
        &[
            (
                roles::SWEEP_BOTTOM,
                &profile,
                Some(Plane::through(Point::ORIGIN, Direction::X)),
            ),
            (
                roles::SWEEP_TOP,
                &profile,
                Some(Plane::through(Point::new(first, second, 0.0), Direction::Y)),
            ),
        ],
    );
    eprintln!(
        "evolved: recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

#[test]
fn a_sweep_surface_traces_each_strip_to_its_edge() {
    // An open L-shaped wire carried by translation up a straight spine: a
    // sheet of one face per wire edge, no caps.
    let mut model = Model::new();
    let profile = path(
        &mut model,
        &[
            Point::new(1.0, 0.0, 0.0),
            Point::ORIGIN,
            Point::new(0.0, 1.0, 0.0),
        ],
    );
    let spine = path(&mut model, &[Point::ORIGIN, Point::new(0.0, 0.0, 2.0)]);
    let built = ogeom_offset::make_sweep_surface(
        &mut model,
        &profile,
        &spine,
        &ogeom_offset::PipeLaw::Fixed,
        T,
    )
    .unwrap();
    let runs = [leg(
        Transform::IDENTITY,
        Point::ORIGIN,
        Point::ORIGIN,
        Vector::new(0.0, 0.0, 1.0),
        (0.0, 2.0),
    )];
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &profile,
        Some(&[&spine]),
        &runs,
        1e-6,
        &[],
    );
    eprintln!(
        "sweep surface: recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

#[test]
fn a_wall_split_by_a_fuse_still_traces_to_its_edge() {
    // A square pipe up the z axis, fused with a box across its +y wall:
    // that wall comes out in two pieces, above and below the box, and the
    // side walls notched. Through the fuse's history each piece is still
    // its edge's.
    let mut model = Model::new();
    let profile = xy_square(&mut model);
    let spine = path(&mut model, &[Point::ORIGIN, Point::new(0.0, 0.0, 3.0)]);
    let built = ogeom_offset::make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &ogeom_offset::PipeLaw::RotationMinimizing,
        ogeom_offset::PipeCorners::Mitre,
        1e-4,
        T,
    )
    .unwrap();
    let block = ogeom_algo::make_box(
        &mut model,
        ogeom_math::Frame::new(Point::new(-1.0, 0.3, 1.0), Direction::Z, Direction::X, T).unwrap(),
        (2.0, 0.7, 1.0),
        T,
    )
    .unwrap()
    .shape;
    let fused = ogeom_bool::fuse(&mut model, &built.shape, &block, T).unwrap();
    let history = built.history.then(&fused.history);
    let runs = [leg(
        Transform::IDENTITY,
        Point::ORIGIN,
        Point::ORIGIN,
        Vector::new(0.0, 0.0, 1.0),
        (0.0, 3.0),
    )];
    let m = check_sweep(
        &model,
        &fused.shape,
        &history,
        &profile,
        None,
        &runs,
        1e-6,
        &[],
    );
    let pieces = explore_unique(&model, &profile, ShapeType::Edge)
        .unwrap()
        .iter()
        .map(|e| {
            history
                .generated(e)
                .iter()
                .filter(|g| model.kind_of(g).unwrap() == ShapeType::Face)
                .count()
        })
        .max()
        .unwrap();
    assert_eq!(pieces, 2, "the split wall's two pieces");
    eprintln!(
        "fused pipe: recorded pieces within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

/// The square wire of half-width `half` about the z axis at height `z`.
fn ring_at(model: &mut Model, half: f64, z: f64) -> Shape {
    let corners = square_corners(
        Point::new(0.0, 0.0, z),
        Vector::new(1.0, 0.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
        half,
    );
    ogeom_algo::make_polygon(model, &corners, true, T)
        .unwrap()
        .shape
}

/// The frustum of squares narrowing from a half-width of one at the base
/// to a half at height two: each point of the base runs straight to its
/// image at the top.
fn frustum_run() -> [Run; 1] {
    [Run {
        motion: Box::new(|s, p| {
            let k = 1.0 - 0.25 * s;
            Point::new(p.x * k, p.y * k, s)
        }),
        range: (0.0, 2.0),
    }]
}

/// Each edge of `later` generates only walls some edge of `first` does.
fn later_section_walls(model: &Model, history: &History, first: &Shape, later: &Shape) {
    let firsts: Vec<Shape> = explore_unique(model, first, ShapeType::Edge)
        .unwrap()
        .iter()
        .flat_map(|e| history.generated(e).to_vec())
        .collect();
    for edge in explore_unique(model, later, ShapeType::Edge).unwrap() {
        let walls: Vec<Shape> = history
            .generated(&edge)
            .iter()
            .filter(|g| model.kind_of(g).unwrap() == ShapeType::Face)
            .cloned()
            .collect();
        assert!(
            !walls.is_empty(),
            "a later section's edge generates no wall"
        );
        for wall in walls {
            assert!(firsts.iter().any(|f| f.is_same(&wall)));
        }
    }
}

#[test]
fn a_ruled_loft_traces_each_wall_to_both_sections() {
    let mut model = Model::new();
    let bottom = ring_at(&mut model, 1.0, 0.0);
    let top = ring_at(&mut model, 0.5, 2.0);
    let built = ogeom_offset::make_loft(&mut model, &bottom, &top, T).unwrap();
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &bottom,
        Some(&[&top]),
        &frustum_run(),
        1e-6,
        &[
            (
                roles::SWEEP_BOTTOM,
                &bottom,
                Some(Plane::through(Point::ORIGIN, Direction::Z)),
            ),
            (
                roles::SWEEP_TOP,
                &top,
                Some(Plane::through(Point::new(0.0, 0.0, 2.0), Direction::Z)),
            ),
        ],
    );
    later_section_walls(&model, &built.history, &bottom, &top);
    eprintln!(
        "ruled loft: recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

#[test]
fn a_skinned_loft_traces_each_wall_to_every_section() {
    // Three squares on one frustum: the skin through them is that frustum.
    let mut model = Model::new();
    let sections = [
        ring_at(&mut model, 1.0, 0.0),
        ring_at(&mut model, 0.75, 1.0),
        ring_at(&mut model, 0.5, 2.0),
    ];
    let tolerance = 1e-4;
    let built = ogeom_offset::make_loft_skinned(&mut model, &sections, tolerance, T).unwrap();
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &sections[0],
        Some(&[&sections[1], &sections[2]]),
        &frustum_run(),
        tolerance * 2.0,
        &[
            (
                roles::SWEEP_BOTTOM,
                &sections[0],
                Some(Plane::through(Point::ORIGIN, Direction::Z)),
            ),
            (
                roles::SWEEP_TOP,
                &sections[2],
                Some(Plane::through(Point::new(0.0, 0.0, 2.0), Direction::Z)),
            ),
        ],
    );
    later_section_walls(&model, &built.history, &sections[0], &sections[1]);
    later_section_walls(&model, &built.history, &sections[0], &sections[2]);
    eprintln!(
        "skinned loft: recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

/// A square pipe up the z axis and on round a quarter arc of radius two
/// turning toward +x, under the rotation-minimizing frame (each leg built
/// exact: a prism, then a revolution about the arc's axis), and its runs.
#[test]
fn a_pipe_down_a_line_and_an_arc_traces_each_leg_to_its_edge() {
    let mut model = Model::new();
    let profile = xy_square(&mut model);
    let up = 2.0;
    use ogeom_geom::Curve3d as _;
    let joint = ogeom_algo::make_vertex(&mut model, Point::new(0.0, 0.0, up)).shape;
    let foot = ogeom_algo::make_vertex(&mut model, Point::ORIGIN).shape;
    let tip = ogeom_algo::make_vertex(&mut model, Point::new(2.0, 0.0, up + 2.0)).shape;
    let line_edge = {
        let curve: ogeom_geom::Curve =
            ogeom_geom::LineCurve::segment(Point::ORIGIN, Point::new(0.0, 0.0, up), T)
                .unwrap()
                .into();
        let range = curve.domain();
        ogeom_algo::make_edge_between(&mut model, curve, range, &foot, &joint, T)
            .unwrap()
            .shape
    };
    let centre = Point::new(2.0, 0.0, up);
    let frame = ogeom_math::Frame::new(centre, Direction::Y, -Direction::X, T).unwrap();
    let circle = ogeom_math::Circle::new(frame, 2.0, T).unwrap();
    let arc = ogeom_algo::make_edge_between(
        &mut model,
        ogeom_geom::CircleCurve::new(circle).into(),
        (0.0, core::f64::consts::FRAC_PI_2),
        &joint,
        &tip,
        T,
    )
    .unwrap()
    .shape;
    let spine = ogeom_algo::make_wire(&mut model, &[line_edge, arc], T)
        .unwrap()
        .shape;
    let built = ogeom_offset::make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &ogeom_offset::PipeLaw::RotationMinimizing,
        ogeom_offset::PipeCorners::Mitre,
        1e-4,
        T,
    )
    .unwrap();
    let runs = [
        leg(
            Transform::IDENTITY,
            Point::ORIGIN,
            Point::ORIGIN,
            Vector::new(0.0, 0.0, 1.0),
            (0.0, up),
        ),
        turn_run(
            Transform::translation(Vector::new(0.0, 0.0, up)),
            Axis::new(centre, Direction::Y),
            core::f64::consts::FRAC_PI_2,
        ),
    ];
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &profile,
        Some(&[&spine]),
        &runs,
        1e-6,
        &[
            (
                roles::SWEEP_BOTTOM,
                &profile,
                Some(Plane::through(Point::ORIGIN, Direction::Z)),
            ),
            (
                roles::SWEEP_TOP,
                &profile,
                Some(Plane::through(Point::new(2.0, 0.0, up + 2.0), Direction::X)),
            ),
        ],
    );
    eprintln!(
        "line and arc pipe: recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

#[test]
fn a_pipe_with_extended_corners_traces_each_run_on_to_its_edge() {
    // The mitred pipe's spine with the corner extended: each leg runs on
    // past the corner until its section clears the other's far side, so
    // its walls lie on the same two prisms run on.
    let mut model = Model::new();
    let profile = xy_square(&mut model);
    let (up, along) = (3.0, 2.0);
    let corner = Point::new(0.0, 0.0, up);
    let spine = path(
        &mut model,
        &[Point::ORIGIN, corner, Point::new(along, 0.0, up)],
    );
    let built = ogeom_offset::make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &ogeom_offset::PipeLaw::RotationMinimizing,
        ogeom_offset::PipeCorners::Extended,
        1e-4,
        T,
    )
    .unwrap();
    let quarter = Transform::rotation(
        Axis::new(Point::ORIGIN, Direction::Y),
        core::f64::consts::FRAC_PI_2,
    );
    let runs = [
        leg(
            Transform::IDENTITY,
            Point::ORIGIN,
            Point::ORIGIN,
            Vector::new(0.0, 0.0, 1.0),
            (0.0, up + 1.0),
        ),
        leg(
            quarter,
            Point::ORIGIN,
            corner,
            Vector::new(1.0, 0.0, 0.0),
            (-1.0, along),
        ),
    ];
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &profile,
        Some(&[&spine]),
        &runs,
        1e-6,
        &[
            (
                roles::SWEEP_BOTTOM,
                &profile,
                Some(Plane::through(Point::ORIGIN, Direction::Z)),
            ),
            (
                roles::SWEEP_TOP,
                &profile,
                Some(Plane::through(Point::new(along, 0.0, up), Direction::X)),
            ),
        ],
    );
    eprintln!(
        "extended pipe: recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

#[test]
fn a_fixed_frame_pipe_traces_each_wall_to_its_edge() {
    // The section carried by translation alone, up the z axis and then
    // up and along x at forty-five degrees: each wall is its edge's prism
    // down either leg.
    let mut model = Model::new();
    let profile = xy_square(&mut model);
    let corner = Point::new(0.0, 0.0, 2.0);
    let end = Point::new(2.0, 0.0, 4.0);
    let spine = path(&mut model, &[Point::ORIGIN, corner, end]);
    let built = ogeom_offset::make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &ogeom_offset::PipeLaw::Fixed,
        ogeom_offset::PipeCorners::Mitre,
        1e-4,
        T,
    )
    .unwrap();
    let slant = Vector::new(1.0, 0.0, 1.0) * core::f64::consts::FRAC_1_SQRT_2;
    let runs = [
        leg(
            Transform::IDENTITY,
            Point::ORIGIN,
            Point::ORIGIN,
            Vector::new(0.0, 0.0, 1.0),
            (0.0, 2.0),
        ),
        leg(
            Transform::IDENTITY,
            Point::ORIGIN,
            corner,
            slant,
            (0.0, corner.distance(end)),
        ),
    ];
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &profile,
        Some(&[&spine]),
        &runs,
        1e-6,
        &[(
            roles::SWEEP_BOTTOM,
            &profile,
            Some(Plane::through(Point::ORIGIN, Direction::Z)),
        )],
    );
    eprintln!(
        "fixed pipe: recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

#[test]
fn a_frenet_pipe_down_a_helix_traces_each_skinned_wall() {
    // A unit square square to a helix of radius ten and pitch five, held
    // in the helix's Frenet frame (its normal pointing at the axis) for
    // one turn: walls skinned to the tolerance.
    use ogeom_geom::Curve3d as _;
    let mut model = Model::new();
    let helix = ogeom_geom::HelixCurve::conical(
        ogeom_math::Frame::WORLD,
        10.0,
        5.0,
        0.0,
        0.0,
        core::f64::consts::TAU,
    )
    .unwrap();
    let range = helix.domain();
    let frame_at = move |t: f64| -> ogeom_math::Frame {
        let at = helix.point_at(t, T).unwrap();
        let tangent = Direction::new(helix.d1_at(t, T).unwrap(), T).unwrap();
        let inward = Direction::new(Vector::new(-at.x, -at.y, 0.0), T).unwrap();
        ogeom_math::Frame::new(at, tangent, inward, T).unwrap()
    };
    let start = frame_at(range.0);
    let corners = square_corners(start.origin(), start.x().vector(), start.y().vector(), 0.5);
    let profile = polygon_face(&mut model, &corners);
    let edge = ogeom_algo::make_edge(&mut model, helix.into(), range, T)
        .unwrap()
        .shape;
    let spine = ogeom_algo::make_wire(&mut model, &[edge], T).unwrap().shape;
    let tolerance = 1e-3;
    let built =
        ogeom_offset::make_pipe_shell(&mut model, &profile, &spine, true, tolerance, T).unwrap();
    let runs = [Run {
        motion: Box::new(move |t, p| {
            (Transform::from_frame(&frame_at(t)) * Transform::to_frame(&start)).apply(p)
        }),
        range,
    }];
    let end = frame_at(range.1);
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &profile,
        Some(&[&spine]),
        &runs,
        // The walls stand a few thousandths off the closed-form Frenet
        // sweep here and the other edges' about one off: the bound tells
        // the strips apart and leaves the skin's accuracy to its own tests.
        1e-2,
        &[
            (
                roles::SWEEP_BOTTOM,
                &profile,
                Some(Plane::through(start.origin(), start.z())),
            ),
            (
                roles::SWEEP_TOP,
                &profile,
                Some(Plane::through(end.origin(), end.z())),
            ),
        ],
    );
    eprintln!(
        "Frenet helix pipe: recorded walls within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}

#[test]
fn an_open_wire_down_a_cornered_spine_traces_each_sheet_strip() {
    // An open L of two edges swept up the z axis and on along x, mitred:
    // a sheet, each strip two pieces, no caps.
    let mut model = Model::new();
    let profile = path(
        &mut model,
        &[
            Point::new(0.5, -0.5, 0.0),
            Point::new(-0.5, -0.5, 0.0),
            Point::new(-0.5, 0.5, 0.0),
        ],
    );
    let (up, along) = (3.0, 2.0);
    let corner = Point::new(0.0, 0.0, up);
    let spine = path(
        &mut model,
        &[Point::ORIGIN, corner, Point::new(along, 0.0, up)],
    );
    let built = ogeom_offset::make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &ogeom_offset::PipeLaw::RotationMinimizing,
        ogeom_offset::PipeCorners::Mitre,
        1e-4,
        T,
    )
    .unwrap();
    let quarter = Transform::rotation(
        Axis::new(Point::ORIGIN, Direction::Y),
        core::f64::consts::FRAC_PI_2,
    );
    let runs = [
        leg(
            Transform::IDENTITY,
            Point::ORIGIN,
            Point::ORIGIN,
            Vector::new(0.0, 0.0, 1.0),
            (0.0, up + 1.0),
        ),
        leg(
            quarter,
            Point::ORIGIN,
            corner,
            Vector::new(1.0, 0.0, 0.0),
            (-1.0, along),
        ),
    ];
    let m = check_sweep(
        &model,
        &built.shape,
        &built.history,
        &profile,
        Some(&[&spine]),
        &runs,
        1e-6,
        &[],
    );
    eprintln!(
        "open wire sheet: recorded strips within {:e}, others at least {:e} off",
        m.worst_on, m.least_off
    );
}
