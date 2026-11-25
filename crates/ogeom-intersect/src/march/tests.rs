use super::split::split_at_touches_as_traced;
use super::*;
use ogeom_geom::{CylinderSurface, PlaneSurface, SphereSurface};
use ogeom_math::{Cylinder, Direction, Frame, Plane, Sphere};

const T: Tolerances = Tolerances::millimetres();

/// A unit circle in the xy plane sampled `per_lap` times a lap for
/// `laps` laps, closed on its start, each sample's chart state its
/// angle.
fn circle_laps(per_lap: usize, laps: usize) -> Traced {
    let n = per_lap * laps;
    let at = |k: usize| {
        #[allow(clippy::cast_precision_loss)]
        let angle = core::f64::consts::TAU * (k % per_lap) as f64 / per_lap as f64;
        (Point::new(angle.cos(), angle.sin(), 0.0), (angle, 0.0))
    };
    let samples: Vec<(Point, (f64, f64))> = (0..=n).map(|k| at(k % n.max(1))).collect();
    Traced {
        points: samples.iter().map(|s| s.0).collect(),
        on_a: samples.iter().map(|s| s.1).collect(),
        on_b: samples.iter().map(|s| s.1).collect(),
        stopped: Stopped::Closed,
    }
}

/// A closed branch traced three times round is held to its first lap,
/// ending on its start; a single lap, and an open branch, are kept as
/// they are.
#[test]
fn a_loop_traced_round_again_is_held_to_its_first_lap() {
    // Each lap's samples shifted a hair along the circle, so the second
    // lap's segment passes over the start rather than through a sample.
    let mut three = circle_laps(64, 3);
    for (k, p) in three
        .points
        .iter_mut()
        .enumerate()
        .skip(64)
        .take(64 * 2 - 1)
    {
        let angle = p.y.atan2(p.x) + 1e-3 * if k < 128 { 1.0 } else { 2.0 };
        *p = Point::new(angle.cos(), angle.sin(), 0.0);
    }
    let held = first_lap(three, 1e-2);
    let lap = circle_laps(64, 1).points.len();
    assert_eq!(held.points.len(), lap, "{} samples", held.points.len());
    assert!(held.points[held.points.len() - 1].distance(held.points[0]) <= 0.0);
    assert_eq!(held.on_a.len(), held.points.len());
    let one = circle_laps(64, 1);
    assert_eq!(first_lap(one.clone(), 1e-2), one);
    let mut open = circle_laps(64, 3);
    open.stopped = Stopped::LeftTheDomain;
    assert_eq!(first_lap(open.clone(), 1e-2), open);
}

fn cylinder(origin: Point, axis: Vector, radius: f64, height: (f64, f64)) -> SurfaceGeometry {
    let frame = Frame::new(
        origin,
        Direction::new(axis, T).unwrap(),
        Direction::from_cross(axis, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
        T,
    )
    .unwrap();
    CylinderSurface::new(Cylinder::new(frame, radius, T).unwrap(), height)
        .unwrap()
        .into()
}

/// A unit drum spelt as its exact clamped patch, which closes round
/// without repeating, cut by a plane touching one of its end rims off
/// the patch's seam: `x cos a + y sin a + c d = 1`, `a` the touch's
/// angle and `d` the depth into the drum from that rim. The section is
/// the ellipse of semi-axes 1 and `b = √(1 + 1/c²)`, touching the rim.
/// A walk from a seed on one side of the touch stops on it, and one
/// from the other side may run over it and round the seam: the arc the
/// first traced is then a stretch of the second, and is returned once,
/// and the stretch from the touch to the seam is returned too. So the
/// arcs' polylines together run as long as the ellipse's perimeter `P`,
/// by Simpson's rule, less what chords sagging at most the march's
/// chord `h` cut off: a chord of a stretch of curvature `κ` sagging `s`
/// is shorter than it by `κ·s·L/3` to leading order, `L` its length, at
/// most `b·h·L/3` here, and at most `b·h·P/3` all round.
#[test]
fn a_plane_touching_a_spline_drums_rim_is_traced_once() {
    let drum: SurfaceGeometry = SurfaceGeometry::from(
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 1.0, T).unwrap(), (0.0, 2.0)).unwrap(),
    )
    .to_bspline(T)
    .unwrap()
    .into();
    for c in [1.5_f64, 2.5, 4.0] {
        let b = (1.0 + 1.0 / (c * c)).sqrt();
        let quarter = core::f64::consts::FRAC_PI_2 / 512.0;
        let perimeter = 4.0 * quarter / 3.0
            * (0..=512)
                .map(|k| {
                    let t = quarter * f64::from(k);
                    let weight = if k == 0 || k == 512 {
                        1.0
                    } else if k % 2 == 1 {
                        4.0
                    } else {
                        2.0
                    };
                    weight * t.sin().hypot(b * t.cos())
                })
                .sum::<f64>();
        for k in 0..12 {
            let a = f64::from(k) / 12.0 * core::f64::consts::TAU + 0.1;
            for (height, up) in [(0.0, 1.0), (2.0, -1.0)] {
                let plane: SurfaceGeometry = PlaneSurface::over(
                    Plane::through(
                        Point::new(a.cos(), a.sin(), height),
                        Direction::new(Vector::new(a.cos(), a.sin(), up * c), T).unwrap(),
                    ),
                    (-8.0, 8.0),
                    (-8.0, 8.0),
                )
                .unwrap()
                .into();
                let found = branches(&drum, &plane, Marching::default(), T).unwrap();
                let length: f64 = found
                    .iter()
                    .flat_map(|branch| branch.points.windows(2))
                    .map(|pair| pair[0].distance(pair[1]))
                    .sum();
                // A polyline through samples of a convex loop runs
                // shorter than the loop. A stretch returned twice adds
                // its length, and one missed takes its length away.
                let cut_off = b * Marching::default().chord * perimeter / 3.0;
                assert!(
                    length < perimeter && perimeter - length <= cut_off,
                    "c {c}, touch at {a} on the rim at {height}: {} arcs run {length}, \
                     the ellipse {perimeter}, chords cut off at most {cut_off}",
                    found.len()
                );
            }
        }
    }
}

/// A drum of radius 0.5 along x, its axis at y = `c` and z = 1.5,
/// ending in the plane x = 2.5 tangent to a bore of radius 2.5 along z.
/// The drum's wall meets the bore in a loop whose x is greatest, 2.5,
/// where it crosses y = 0: at z = 1.5 ± √(0.25 − c²), two points of the
/// drum's end rim, which the loop touches and turns back from. So the
/// loop is two arcs between the touches, each ending on both, with
/// either surface first and the rim at either end of the drum's height.
/// At c = 0.45 and 0.475 the arc below y = 0 stands within 0.0005 and
/// 0.00013 of the rim. A walk ending short of a touch or running over
/// it leaves an arc's end microns to millimetres off the rim's point,
/// or the loop whole.
#[test]
fn a_section_touching_a_drum_rim_ends_on_the_touch() {
    let bore: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 2.5, T).unwrap(), (0.0, 3.0))
            .unwrap()
            .into();
    for c in [0.2, 0.45, 0.475] {
        let lift = (0.25_f64 - c * c).sqrt();
        let touches = [
            Point::new(2.5, 0.0, 1.5 - lift),
            Point::new(2.5, 0.0, 1.5 + lift),
        ];
        for (origin, axis) in [
            (Point::new(1.0, c, 1.5), Direction::X),
            (
                Point::new(2.5, c, 1.5),
                Direction::new(-Vector::X, T).unwrap(),
            ),
        ] {
            let frame = Frame::new(origin, axis, Direction::Z, T).unwrap();
            let drum: SurfaceGeometry =
                CylinderSurface::new(Cylinder::new(frame, 0.5, T).unwrap(), (0.0, 1.5))
                    .unwrap()
                    .into();
            for (a, b) in [(&bore, &drum), (&drum, &bore)] {
                let found = branches(a, b, Marching::default(), T).unwrap();
                assert_eq!(found.len(), 2, "{c}, {axis:?}: the loop's two arcs");
                for touch in touches {
                    let ends = found
                        .iter()
                        .flat_map(|arc| [arc.points[0], arc.points[arc.points.len() - 1]])
                        .filter(|end| end.distance(touch) <= T.confusion())
                        .count();
                    assert_eq!(ends, 2, "{c}, {axis:?}: the arcs ending on {touch:?}");
                }
            }
        }
    }
}

/// Two drums across each other, the second's foot a little above the
/// first's axis and the first's end a little past the second's axis.
/// Their crossing is one loop, and the domains cut it into two arcs:
/// a long one, and one 0.09 long between the first drum's rim and the
/// second's, inside one pair of cells and crossing neither. It is
/// seeded where the rims pierce the other drum.
#[test]
fn a_short_arc_cut_off_between_two_rims_is_found() {
    let drum = |origin: Point, axis: Direction, x: Direction, radius: f64, height: f64| {
        let frame = Frame::new(origin, axis, x, T).unwrap();
        SurfaceGeometry::from(
            CylinderSurface::new(Cylinder::new(frame, radius, T).unwrap(), (0.0, height)).unwrap(),
        )
    };
    let a = drum(
        Point::new(-1.419_880_124_370_016_3, -1.0, 3.0),
        Direction::X,
        Direction::Y,
        1.773_938_116_486_465_8,
        2.5,
    );
    let b = drum(
        Point::new(0.510_316_713_993_857_4, 1.5, 2.230_205_316_817_85),
        Direction::Z,
        Direction::X,
        1.040_566_819_682_834_4,
        4.0,
    );
    let found = branches(&a, &b, Marching::default(), T).unwrap();
    // The short arc runs between the two rims, all of it past
    // x = 1.02; the long one reaches back to x = -0.24.
    assert!(
        found
            .iter()
            .any(|branch| branch.points.iter().all(|p| p.x > 1.02)),
        "{} branches, none the short arc",
        found.len()
    );
}

/// The patch over `0 ≤ x ≤ 1`, `0 ≤ y ≤ 2` standing `z = s y + q y²`
/// off the plane `z = 0`: linear along `x`, its quadratic in `y` with
/// the middle control point on the start's heading half way. Its
/// border `y = 0` lies in the plane, and it leaves the plane there at
/// the sine `s / √(1 + s²)`, bowing on at `q`.
fn leaving(s: f64, q: f64) -> SurfaceGeometry {
    let points: Vec<Point> = [0.0, 1.0]
        .into_iter()
        .flat_map(|x| {
            [(0.0, 0.0), (1.0, s), (2.0, 2.0 * s + 4.0 * q)].map(|(y, z)| Point::new(x, y, z))
        })
        .collect();
    ogeom_geom::BSplineSurface::new(
        ogeom_math::KnotVector::clamped_uniform(1, 2).unwrap(),
        ogeom_math::KnotVector::clamped_uniform(2, 3).unwrap(),
        &ogeom_math::ControlGrid::new(points, 2, 3).unwrap(),
        T,
    )
    .unwrap()
    .into()
}

/// A border lying in the other surface is a stretch of the section
/// however shallowly the surface leaves it, so long as the angle takes
/// it off: a sine `s` puts it `s d` off a distance `d` in, past the
/// confusion distance by the grid's first row in, here `d = 1/12`. One
/// leaving tangentially, a sine of naught there and rising only as it
/// bows, is a contact and gives no stretch, and so is one leaving at
/// `2e-4` (the little a fit of a tangent blend leaves) while bowing at a
/// curvature of one, which stands it off two hundred times further there.
#[test]
fn a_border_left_shallowly_is_a_stretch_and_one_left_tangentially_is_not() {
    let ground = plane(Point::ORIGIN, Vector::Z);
    let options = Marching::default();
    for (s, q, stretches) in [
        (0.002, 0.005, 1),
        (0.01, 0.0, 1),
        (0.02, 0.02, 1),
        (0.0, 0.02, 0),
        (2e-4, 0.5, 0),
    ] {
        let patch = leaving(s, q);
        for (a, b) in [(&ground, &patch), (&patch, &ground)] {
            let mesh = (sample(a, options.grid, T), sample(b, options.grid, T));
            let found = border_stretches(a, b, (&mesh.0, &mesh.1), options, T).unwrap();
            assert_eq!(found.len(), stretches, "sine {s}, bow {q}: {found:?}");
            for piece in &found {
                let ends = [piece.points[0], piece.points[piece.points.len() - 1]];
                for p in &piece.points {
                    assert!(p.y.abs().max(p.z.abs()) <= T.confusion(), "sine {s}: {p:?}");
                }
                let mut xs = ends.map(|p| p.x);
                xs.sort_by(f64::total_cmp);
                assert!(
                    xs[0].abs() <= T.confusion() && (xs[1] - 1.0).abs() <= T.confusion(),
                    "sine {s}: the stretch runs {xs:?}"
                );
            }
        }
    }
}

fn sphere(centre: Point, radius: f64) -> SurfaceGeometry {
    SphereSurface::new(Sphere::centred(centre, radius, T).unwrap()).into()
}

fn plane(origin: Point, normal: Vector) -> SurfaceGeometry {
    PlaneSurface::over(
        Plane::through(origin, Direction::new(normal, T).unwrap()),
        (-8.0, 8.0),
        (-8.0, 8.0),
    )
    .unwrap()
    .into()
}

/// How far a point is from a quadric, in closed form.
fn off(surface: &SurfaceGeometry, p: Point) -> f64 {
    match surface {
        SurfaceGeometry::Plane(x) => x.plane().distance_to(p),
        SurfaceGeometry::Sphere(x) => x.sphere().distance_to(p),
        SurfaceGeometry::Cylinder(x) => x.cylinder().distance_to(p),
        _ => 0.0,
    }
}

/// The worst distance from a traced branch to either surface.
fn deviation(a: &SurfaceGeometry, b: &SurfaceGeometry, traced: &Traced) -> f64 {
    traced
        .points
        .iter()
        .map(|p| off(a, *p).abs().max(off(b, *p).abs()))
        .fold(0.0_f64, f64::max)
}

#[test]
fn a_plane_through_a_bent_strip_seeds_both_branches() {
    // A cubic strip bent into an arch crosses a level plane twice. The
    // plane's domain is the unbounded one a face carries, and its
    // sampling cells span a million units; merging seeds at *its*
    // spacing would call the two crossings one branch and trace only
    // one of them. Seeds merge at the finer surface's spacing instead.
    use ogeom_geom::BSplineSurface;
    use ogeom_math::{ControlGrid, KnotVector};
    let mut points = Vec::new();
    for i in 0..7 {
        let a = core::f64::consts::PI * f64::from(i) / 6.0;
        for j in 0..2 {
            points.push(Point::new(2.0 * a.cos(), f64::from(j), 2.0 * a.sin()));
        }
    }
    let grid = ControlGrid::new(points, 7, 2).unwrap();
    let strip: SurfaceGeometry = BSplineSurface::new(
        KnotVector::clamped_uniform(3, 7).unwrap(),
        KnotVector::clamped_uniform(1, 2).unwrap(),
        &grid,
        T,
    )
    .unwrap()
    .into();
    let level: SurfaceGeometry = PlaneSurface::over(
        Plane::through(Point::new(0.0, 0.0, 1.0), Direction::Z),
        (-1.0e9, 1.0e9),
        (-1.0e9, 1.0e9),
    )
    .unwrap()
    .into();
    let options = Marching {
        chord: 1e-5,
        ..Marching::default()
    };
    let found = branches(&strip, &level, options, T).unwrap();
    assert_eq!(
        found.len(),
        2,
        "the arch crosses the level twice: {}",
        found.len()
    );
    for branch in &found {
        assert!(!branch.closed());
        for p in &branch.points {
            assert!((p.z - 1.0).abs() < 1e-4, "on the level: {p:?}");
        }
    }
}

#[test]
fn two_crossed_cylinders_are_traced_onto_both_of_them() {
    // The case with no closed form, and the one the analytic module
    // explicitly refuses: two cylinders on perpendicular axes meet in a
    // quartic space curve. Every point of the trace must be on both.
    let a = cylinder(Point::ORIGIN, Vector::Z, 1.0, (-4.0, 4.0));
    let b = cylinder(Point::ORIGIN, Vector::X, 1.0, (-4.0, 4.0));
    let options = Marching {
        chord: 1e-5,
        ..Marching::default()
    };

    // Two equal cylinders crossing at right angles meet in two ellipses,
    // the Steinmetz solid's seams, crossing where the cylinders touch at
    // (0, -1, 0) and (0, 1, 0): four arcs, each from one touch to the
    // other.
    let found = branches(&a, &b, options, T).unwrap();
    assert_eq!(found.len(), 4, "the two ellipses' halves");

    let touches = [Point::new(0.0, -1.0, 0.0), Point::new(0.0, 1.0, 0.0)];
    let mut worst = 0.0_f64;
    for branch in &found {
        let ends = [branch.points[0], branch.points[branch.points.len() - 1]];
        for touch in touches {
            assert!(
                ends.iter().any(|end| end.distance(touch) < 1e-9),
                "an arc from {:?} to {:?} misses the touch {touch:?}",
                ends[0],
                ends[1]
            );
        }
        assert!(
            branch.points.len() > 100,
            "a branch of only {} points",
            branch.points.len()
        );
        worst = worst.max(deviation(&a, &b, branch));
    }
    println!(
        "crossed cylinders: {} branches, worst deviation {worst:e}",
        found.len()
    );
    assert!(worst < 1e-7, "traced off the surfaces by {worst:e}");
}

/// A cylinder lying inside a torus's outer equator and touching it at
/// (13, 0, 0), where the torus's chart has its corner, meets it in a
/// figure eight. Traced either way round, the section comes back as the
/// two loops' four halves, each ending on the touch, and the loops'
/// length is the section's.
#[test]
fn a_figure_eight_is_cut_at_its_double_point() {
    let radius = 4.353_623_591_855_474;
    let torus: SurfaceGeometry =
        ogeom_geom::TorusSurface::new(ogeom_math::Torus::new(Frame::WORLD, 10.0, 3.0, T).unwrap())
            .into();
    let frame = Frame::new(
        Point::new(13.0 - radius, 0.0, -10.0),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let drill: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(frame, radius, T).unwrap(), (0.0, 20.0))
            .unwrap()
            .into();
    let options = Marching {
        chord: 6e-6,
        ..Marching::default()
    };
    // The section's length: the drill's wall unrolled, z the root of 9
    // less the square of the distance from the tube's centre circle, at
    // each angle round the drill where that is positive, four times
    // over for the two loops' upper and lower halves.
    let steps = 200_000;
    let mut expected = 0.0;
    let at = |t: f64| {
        let (x, y) = (13.0 - radius + radius * t.cos(), radius * t.sin());
        let off = x.hypot(y) - 10.0;
        ((9.0 - off * off).max(0.0).sqrt(), x, y)
    };
    for k in 0..steps {
        let (t0, t1) = (
            core::f64::consts::TAU * f64::from(k) / f64::from(steps),
            core::f64::consts::TAU * f64::from(k + 1) / f64::from(steps),
        );
        let ((z0, x0, y0), (z1, x1, y1)) = (at(t0), at(t1));
        if z0 > 0.0 || z1 > 0.0 {
            expected += 2.0 * Point::new(x0, y0, z0).distance(Point::new(x1, y1, z1));
        }
    }
    let touch = Point::new(13.0, 0.0, 0.0);
    for (a, b) in [(&torus, &drill), (&drill, &torus)] {
        let found = branches(a, b, options, T).unwrap();
        assert_eq!(found.len(), 4, "the figure eight's four halves");
        let mut length = 0.0;
        for branch in &found {
            let ends = [branch.points[0], branch.points[branch.points.len() - 1]];
            assert!(
                ends.iter().filter(|end| end.distance(touch) < 1e-9).count() == 1,
                "a half from {:?} to {:?} does not end once on the touch",
                ends[0],
                ends[1]
            );
            length += branch
                .points
                .windows(2)
                .map(|w| w[0].distance(w[1]))
                .sum::<f64>();
            assert!(deviation(a, b, branch) < 1e-6, "traced off the surfaces");
        }
        assert!(
            (length - expected).abs() < 1e-3,
            "the halves run {length} against the section's {expected}"
        );
    }
}

#[test]
fn unequal_crossed_cylinders_meet_in_two_curves_as_well() {
    // The non-degenerate cousin. Equal radii put the two curves through
    // each other at the tangency points, so getting the count right there
    // says less than getting it right here, where there is no singularity
    // for a tracer to be lucky about.
    let a = cylinder(Point::ORIGIN, Vector::Z, 1.0, (-4.0, 4.0));
    let b = cylinder(Point::ORIGIN, Vector::X, 1.6, (-4.0, 4.0));
    let options = Marching {
        chord: 1e-5,
        ..Marching::default()
    };

    let found = branches(&a, &b, options, T).unwrap();
    assert_eq!(found.len(), 2);
    for branch in &found {
        assert!(branch.closed());
        assert!(deviation(&a, &b, branch) < 1e-7);
    }
}

#[test]
fn a_traced_circle_agrees_with_the_circle_it_should_be() {
    // A sphere cut by a plane through its centre is a circle of known
    // radius, and the marcher does not know that. Tracing it and checking
    // against the closed form is the strongest single check there is: it
    // tests the tracer against an answer derived independently of it.
    let s = sphere(Point::ORIGIN, 3.0);
    let cut = plane(Point::ORIGIN, Vector::Z);
    let options = Marching {
        chord: 1e-6,
        ..Marching::default()
    };

    let found = seeds(&s, &cut, options, T).unwrap();
    assert!(!found.is_empty());
    let branch = trace(&s, &cut, found[0], options, T).unwrap();

    assert!(
        branch.closed(),
        "a plane through a sphere gives a closed loop"
    );
    for p in &branch.points {
        let radius = (p.x * p.x + p.y * p.y).sqrt();
        assert!(
            (radius - 3.0).abs() < 1e-7,
            "a point at radius {radius} on a circle of 3"
        );
        assert!(p.z.abs() < 1e-7, "off the cutting plane by {}", p.z);
    }
}

/// A branch the backward half closes is the loop, whatever stopped the
/// forward half. An oblique plane cuts a drum in an ellipse, which takes
/// a different number of steps round each way from most seeds. From one
/// that goes round in fewer backwards, with points enough only for that,
/// the forward half runs out and the backward one closes: joined, the
/// two ran past the seed a second time and were reported truncated.
#[test]
fn a_branch_closed_going_backwards_is_the_loop() {
    let drum = cylinder(Point::new(0.0, 0.0, -4.0), Vector::Z, 1.0, (0.0, 8.0));
    let cut = plane(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0));
    let options = Marching {
        chord: 1e-5,
        ..Marching::default()
    };
    let found = seeds(&drum, &cut, options, T).unwrap();
    let round = trace(&drum, &cut, found[0], options, T).unwrap();
    let Some((from, ahead, behind)) = (0..round.points.len())
        .step_by(37)
        .map(|k| {
            let from = Contact {
                on_a: round.on_a[k],
                on_b: round.on_b[k],
                point: round.points[k],
            };
            let ahead = walk(&drum, &cut, from, 1.0, options, T).unwrap();
            let behind = walk(&drum, &cut, from, -1.0, options, T).unwrap();
            (from, ahead, behind)
        })
        .find(|(_, ahead, behind)| {
            ahead.stopped == Stopped::Closed
                && behind.stopped == Stopped::Closed
                && behind.points.len() < ahead.points.len()
        })
    else {
        panic!("a seed the backward walk closes from in fewer points");
    };

    let options = Marching {
        max_points: behind.points.len(),
        ..options
    };
    let branch = trace(&drum, &cut, from, options, T).unwrap();
    assert_eq!(branch.stopped, Stopped::Closed, "the backward half closed");
    assert_eq!(branch.points.len(), behind.points.len(), "once round");
    assert_eq!(branch.points[0], branch.points[branch.points.len() - 1]);
    assert!(
        (branch.points[1] - from.point).dot(ahead.points[1] - from.point) > 0.0,
        "it runs the way the forward half set off"
    );
}

#[test]
fn a_branch_that_leaves_the_surface_says_so_rather_than_stopping_quietly() {
    // A truncated branch and a finished one are different answers, and a
    // caller that cannot tell them apart treats the first as the second.
    let s = sphere(Point::ORIGIN, 3.0);
    let cut = plane(Point::new(0.0, 0.0, 0.0), Vector::Z);
    let options = Marching {
        chord: 1e-4,
        max_points: 8,
        ..Marching::default()
    };
    let found = seeds(&s, &cut, options, T).unwrap();
    let branch = trace(&s, &cut, found[0], options, T).unwrap();
    assert_eq!(branch.stopped, Stopped::RanOut);
    assert!(!branch.complete(), "a truncated branch is not complete");
}

#[test]
fn tangent_surfaces_are_refused_rather_than_followed_onto_a_guess() {
    // Where the normals are parallel the intersection has no single
    // direction, and marching through such a point is how a tracer changes
    // branch without noticing. A sphere resting on a plane is the case.
    let s = sphere(Point::new(0.0, 0.0, 3.0), 3.0);
    let ground = plane(Point::ORIGIN, Vector::Z);
    let touch = Contact {
        on_a: (0.0, -core::f64::consts::FRAC_PI_2),
        on_b: (0.0, 0.0),
        point: Point::ORIGIN,
    };
    let err = trace(&s, &ground, touch, Marching::default(), T).unwrap_err();
    assert!(err.to_string().contains("tangent"), "unexpected: {err}");
}

#[test]
fn the_number_of_branches_is_the_number_there_are() {
    // The failure the accuracy measure cannot see. Every point of one
    // circle is on both surfaces, so returning one of two scores perfectly;
    // the count is the only thing that catches it.
    let options = Marching {
        chord: 1e-5,
        ..Marching::default()
    };

    // A sphere cut by a plane off its centre: one circle.
    let one = branches(
        &sphere(Point::ORIGIN, 3.0),
        &plane(Point::new(0.0, 0.0, 1.0), Vector::Z),
        options,
        T,
    )
    .unwrap();
    assert_eq!(one.len(), 1, "one plane through a sphere cuts one circle");
    assert!(one[0].closed());

    // A sphere and a coaxial cylinder narrower than it: two circles, one
    // above and one below.
    let two = branches(
        &sphere(Point::ORIGIN, 3.0),
        &cylinder(Point::ORIGIN, Vector::Z, 1.5, (-4.0, 4.0)),
        options,
        T,
    )
    .unwrap();
    assert_eq!(two.len(), 2, "a coaxial cylinder cuts a sphere twice");
    for branch in &two {
        assert!(branch.closed(), "each is a closed circle");
    }
    // And they are on opposite sides, rather than the same one twice.
    let heights: Vec<f64> = two.iter().map(|b| b.points[0].z).collect();
    assert!(
        heights[0] * heights[1] < 0.0,
        "both branches came back on the same side: {heights:?}"
    );
}

/// A thin drum at right angles to a long drill, its axis passing 46
/// from the drill's (radius 56), runs in through the drill's wall and
/// out again: two loops, about the two points where the drum's axis
/// pierces the wall, 64 apart. That is closer than a drill cell's long
/// side, and nearer than either surface's length over the grid.
#[test]
fn a_thin_drum_through_a_long_drill_s_wall_cuts_it_twice() {
    let (radius, offset) = (56.0, 46.0);
    let drill = cylinder(
        Point::new(-900.0, 0.0, 0.0),
        Vector::X,
        radius,
        (0.0, 1800.0),
    );
    let drum = cylinder(
        Point::new(5.0, 0.0, offset),
        Vector::Y,
        5.0,
        (-1000.0, 1000.0),
    );
    let found = branches(&drum, &drill, Marching::default(), T).unwrap();
    assert_eq!(found.len(), 2, "the drum passes through the wall twice");
    let pierce = (radius * radius - offset * offset).sqrt();
    let mut sides: Vec<f64> = found
        .iter()
        .map(|branch| {
            assert!(branch.closed(), "each is a loop round the drum");
            let mean = branch.points.iter().map(|p| p.y).sum::<f64>();
            #[allow(clippy::cast_precision_loss)]
            let mean = mean / branch.points.len() as f64;
            assert!(
                (mean.abs() - pierce).abs() < 5.0,
                "a loop about y = {mean}, the axis pierces at {pierce}"
            );
            mean
        })
        .collect();
    sides.sort_by(f64::total_cmp);
    assert!(sides[0] < 0.0 && sides[1] > 0.0, "{sides:?}");
}

#[test]
fn a_branch_thinner_than_the_sampling_is_found_from_the_surfaces_extremum() {
    // Two spheres barely overlapping meet in a circle far smaller than a
    // coarse grid's cell: no pair of cells crosses on it, and a sphere
    // has no border to pierce. The approach of the two surfaces inside
    // the circle seeds it at any grid, and a finer grid still finds it
    // by crossing cells.
    let a = sphere(Point::ORIGIN, 3.0);
    let b = sphere(Point::new(5.98, 0.0, 0.0), 3.0);
    let radius = (9.0_f64 - 2.99 * 2.99).sqrt();
    for grid in [6, 24, 120] {
        let options = Marching {
            grid,
            ..Marching::default()
        };
        assert!(
            !seeds(&a, &b, options, T).unwrap().is_empty(),
            "grid {grid}"
        );
        let found = branches(&a, &b, options, T).unwrap();
        assert_eq!(found.len(), 1, "one circle at grid {grid}");
        assert!(found[0].closed());
        assert!(deviation(&a, &b, &found[0]) < 1e-7);
        for p in &found[0].points {
            assert!((p.x - 2.99).abs() < 1e-6, "in the mid plane: {p:?}");
            assert!(
                (p.y.hypot(p.z) - radius).abs() < 1e-6,
                "on the circle: {p:?}"
            );
        }
    }
}

/// A bicubic patch over the square `[0, 10]`, level at z = 0 but for a
/// bump raised by one control point near the middle, its top at
/// `peak`.
fn bumped(peak: f64) -> SurfaceGeometry {
    use ogeom_geom::BSplineSurface;
    use ogeom_math::{ControlGrid, KnotVector};
    const N: usize = 41;
    // A raised control point lifts its patch by four ninths of its own
    // height: the cubic basis peaks at two thirds each way.
    let lift = peak * 9.0 / 4.0;
    let mut points = Vec::new();
    for i in 0..N {
        for j in 0..N {
            #[allow(clippy::cast_precision_loss)]
            let (x, y) = (10.0 * i as f64 / 40.0, 10.0 * j as f64 / 40.0);
            let z = if (i, j) == (19, 22) { lift } else { 0.0 };
            points.push(Point::new(x, y, z));
        }
    }
    BSplineSurface::new(
        KnotVector::clamped_uniform(3, N).unwrap(),
        KnotVector::clamped_uniform(3, N).unwrap(),
        &ControlGrid::new(points, N, N).unwrap(),
        T,
    )
    .unwrap()
    .into()
}

/// A bilinear patch level at `height`, wider than the bumped patch.
fn level_patch(height: f64) -> SurfaceGeometry {
    use ogeom_geom::BSplineSurface;
    use ogeom_math::{ControlGrid, KnotVector};
    let points = vec![
        Point::new(-1.0, -1.0, height),
        Point::new(-1.0, 11.0, height),
        Point::new(11.0, -1.0, height),
        Point::new(11.0, 11.0, height),
    ];
    BSplineSurface::new(
        KnotVector::clamped_uniform(1, 2).unwrap(),
        KnotVector::clamped_uniform(1, 2).unwrap(),
        &ControlGrid::new(points, 2, 2).unwrap(),
        T,
    )
    .unwrap()
    .into()
}

#[test]
fn a_small_bump_through_a_spline_patch_cuts_a_loop_inside_both() {
    // The bump's top pokes a hundredth of its height through a level
    // patch: a closed loop far smaller than a grid cell, touching no
    // border of either patch.
    let peak = 0.5;
    let bump = bumped(peak);
    let level = level_patch(peak * 0.99);
    let found = branches(&bump, &level, Marching::default(), T).unwrap();
    assert_eq!(found.len(), 1, "one loop round the bump's top");
    let loop_ = &found[0];
    assert!(loop_.closed());
    let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
    for (p, &(u, v)) in loop_.points.iter().zip(&loop_.on_a) {
        assert!((p.z - peak * 0.99).abs() < 1e-7, "on the level: {p:?}");
        assert!(bump.point_at(u, v, T).unwrap().distance(*p) < 1e-7);
        low = low.min(p.x);
        high = high.max(p.x);
    }
    assert!(
        high - low < 10.0 / 24.0,
        "thinner than a cell: {}",
        high - low
    );

    // Just touching: the top meets the level within tolerance, and a
    // touch has no curve to follow.
    let touching = level_patch(peak);
    assert!(
        branches(&bump, &touching, Marching::default(), T)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn a_drill_grazing_a_ball_cuts_one_small_loop() {
    // A drill along x passes a thousandth deeper than tangent to a ball:
    // a loop a few hundredths across, inside one cell of either surface.
    let ball = sphere(Point::ORIGIN, 3.0);
    let drill = cylinder(Point::new(0.0, 0.0, 3.999), Vector::X, 1.0, (-4.0, 4.0));
    let found = branches(&drill, &ball, Marching::default(), T).unwrap();
    assert_eq!(found.len(), 1, "one loop where the drill dips in");
    assert!(found[0].closed());
    assert!(deviation(&drill, &ball, &found[0]) < 1e-7);
    for p in &found[0].points {
        assert!(p.x.abs() < 0.2 && p.y.abs() < 0.1, "near the graze: {p:?}");
    }
}

#[test]
fn settings_that_could_not_work_are_refused() {
    let a = sphere(Point::ORIGIN, 1.0);
    let b = plane(Point::ORIGIN, Vector::Z);
    for options in [
        Marching {
            chord: 0.0,
            ..Marching::default()
        },
        Marching {
            grid: 1,
            ..Marching::default()
        },
        Marching {
            max_points: 1,
            ..Marching::default()
        },
    ] {
        assert!(seeds(&a, &b, options, T).is_err());
    }
}

/// The rim of a cone `ρ = 3.2 − z / 10` about the `z` axis rounded at
/// radius ½, over `turn` from where it touches the cone, swept along
/// `y` over `sweep`, turned `angle` about the cone's axis, and the cone:
/// the band stands off the cone by about `0.164 y² − u² / 4` near the
/// node at `u = y = 0`, where the section crosses itself. With the
/// band's centre `c` and its point at `(u, y)` unturned.
#[allow(clippy::type_complexity)]
fn band_and_cone(
    turn: f64,
    sweep: (f64, f64),
    angle: f64,
) -> (
    SurfaceGeometry,
    SurfaceGeometry,
    Box<dyn Fn(f64, f64) -> Point>,
) {
    let lean = 1.01_f64.sqrt();
    let c = lean.mul_add(-0.5, 3.05);
    let frame = Frame::from_axes(
        Point::new(c, 0.0, 1.5),
        Direction::new(Vector::new(1.0, 0.0, 0.1), T).unwrap(),
        Direction::new(Vector::new(-0.1, 0.0, 1.0), T).unwrap(),
        -Direction::Y,
        T,
    )
    .unwrap();
    let arc: ogeom_geom::Curve = ogeom_geom::TrimmedCurve::new(
        ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(frame, 0.5, T).unwrap()).into(),
        0.0,
        turn,
        T,
    )
    .unwrap()
    .into();
    let band: SurfaceGeometry = ogeom_geom::ExtrusionSurface::over(arc, Direction::Y, sweep)
        .unwrap()
        .into();
    let about = ogeom_math::Transform::rotation(
        ogeom_math::Axis {
            location: Point::ORIGIN,
            direction: Direction::Z,
        },
        angle,
    );
    let band = ogeom_geom::Transformable::transformed(&band, &about, T).unwrap();
    let frame = Frame::new(Point::new(0.0, 0.0, 2.0), -Direction::Z, -Direction::X, T).unwrap();
    let cone = ogeom_geom::ConeSurface::new(
        ogeom_math::Cone::new(frame, 3.0, 0.1_f64.atan(), T).unwrap(),
        (0.0, 2.0),
    )
    .unwrap()
    .into();
    let at = move |u: f64, y: f64| {
        let (cos, sin) = (u.cos(), u.sin());
        Point::new(
            0.5f64.mul_add((cos - 0.1 * sin) / lean, c),
            y,
            0.5f64.mul_add((0.1f64.mul_add(cos, sin)) / lean, 1.5),
        )
    };
    (band, cone, Box::new(at))
}

/// Branches cut where touches crowd, as traced: the band and cone of
/// [`band_and_cone`] with the branch to the lower end of the sweep
/// 0.0008 long, under ten chords, in both orders, turned every fifteen
/// degrees. The two branches traced from the seeds are joined at the
/// node as a walk turning the corner there, and cut as traced. Each arc
/// leaving the node must be certified a crossing there, every sample
/// must lie on both surfaces, and an arc must run from the node to
/// where the lower branch crosses the end of the sweep, at the turn `u`
/// solving `√(x(u)² + y²) = 3.2 − z(u) / 10`. The branch was not cut,
/// the surfaces crossing too shallowly along it to count the node; cut,
/// the arm went as debris under ten chords; and no arc was certified.
#[test]
fn crowded_touches_cut_as_traced_keep_a_short_arm_as_a_crossing() {
    let low = -0.000_8;
    let off = |p: Point| p.x.hypot(p.y) - 0.1f64.mul_add(-p.z, 3.2);
    let options = Marching::default();
    for k in 0..24 {
        let angle = (f64::from(k) * 15.0).to_radians();
        let (band, cone, at) = band_and_cone(0.425, (low, 0.02), angle);
        let crossing = {
            let (mut lo, mut hi) = (1e-9, 0.425);
            for _ in 0..200 {
                let mid = f64::midpoint(lo, hi);
                if off(at(mid, low)).signum() == off(at(lo, low)).signum() {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            band.point_at(f64::midpoint(lo, hi), low, T).unwrap()
        };
        let node = band.point_at(0.0, 0.0, T).unwrap();
        for (a, b) in [(&band, &cone), (&cone, &band)] {
            let seeds = seeded(a, b, options, T).unwrap().found;
            let traced: Vec<Traced> = seeds
                .into_iter()
                .filter_map(|seed| trace(a, b, seed, options, T).ok())
                .collect();
            // The two branches, each run from the node out, joined
            // there as a walk turning the corner onto the other arm.
            let from_node = |towards: Point| {
                let mut branch = traced
                    .iter()
                    .find(|t| {
                        let m = t.points.len();
                        [t.points[0], t.points[m - 1]]
                            .iter()
                            .any(|p| p.distance(towards) <= options.chord)
                    })
                    .unwrap_or_else(|| panic!("turned {k} × 15°: nothing traced"))
                    .clone();
                if branch.points[0].distance(node) > branch.points[1].distance(node) {
                    branch = reversed(branch);
                }
                branch
            };
            let top = traced
                .iter()
                .flat_map(|t| [t.points[0], t.points[t.points.len() - 1]])
                .max_by(|p, q| p.distance(node).total_cmp(&q.distance(node)))
                .unwrap_or_else(|| panic!("turned {k} × 15°: nothing traced"));
            let upper = from_node(top);
            let lower = from_node(crossing);
            let mut corner = reversed(upper);
            corner.points.extend(&lower.points[1..]);
            corner.on_a.extend(&lower.on_a[1..]);
            corner.on_b.extend(&lower.on_b[1..]);
            corner.stopped = Stopped::LeftTheDomain;
            let arcs = split_at_touches_as_traced(vec![corner], a, b, options, T);
            let mut reaches = false;
            for (arc, crossing_there) in &arcs {
                for p in &arc.points {
                    let unturned = Point::new(
                        p.x.mul_add(angle.cos(), p.y * angle.sin()),
                        p.y.mul_add(angle.cos(), -p.x * angle.sin()),
                        p.z,
                    );
                    let lean = 1.01_f64.sqrt();
                    let from_band =
                        ((unturned.x - lean.mul_add(-0.5, 3.05)).hypot(unturned.z - 1.5) - 0.5)
                            .abs();
                    assert!(
                        from_band <= T.confusion() && off(*p).abs() <= T.confusion(),
                        "turned {k} × 15°: {p:?} is off the band or the cone"
                    );
                }
                let m = arc.points.len();
                let ends = [arc.points[0], arc.points[m - 1]];
                if ends.iter().any(|p| p.distance(node) <= options.chord) {
                    assert!(crossing_there, "turned {k} × 15°: {arc:?}");
                    reaches |= ends.iter().any(|p| p.distance(crossing) <= options.chord);
                }
            }
            assert!(
                reaches,
                "turned {k} × 15°: no arc runs from the node to {crossing:?}: {arcs:?}"
            );
        }
    }
}
