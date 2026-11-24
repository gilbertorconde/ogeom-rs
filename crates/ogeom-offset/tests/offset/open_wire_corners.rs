//! Open wires swept down spines of lines and arcs under the
//! rotation-minimizing frame and the laws carrying the section as it does,
//! round mitred and rounded corners, against closed forms and against the
//! walls of the closed wire or face they halve.
//!
//! A profile point at `α` outside a corner's axis (along the side the
//! corner turns away from), turning a corner through `θ`: mitred, its path
//! runs on `α tan(θ/2)` along each leg (negative inside the turn); rounded,
//! an inner point's path stops short by as much on each leg and an outer
//! point's turns an arc `α θ` about the axis. A sheet's area is the
//! profile's length times the legs' and those gains integrated along the
//! profile; its boundary is the profile twice and the paths of the
//! profile's two ends.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{FRAC_PI_2, FRAC_PI_3, PI};

use ogeom_algo::{check, curve_length, make_polygon, make_prism, make_wire, surface_properties};
use ogeom_core::{FastMap, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Axis, Direction, Point, Transform, Vector};
use ogeom_mesh::Deflection;
use ogeom_offset::{PipeCorners, PipeLaw, make_pipe_shell_with};
use ogeom_topo::{EdgeRepr, Filter, Location, Model, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

/// The surfaces are planes, cylinders and the planes and cones of a turned
/// line, measured exactly; the bound covers rounding in the trims' cut
/// points.
const EXACT: f64 = 1e-9;

fn sweep(model: &mut Model, profile: &Shape, spine: &Shape, corners: PipeCorners) -> Shape {
    make_pipe_shell_with(
        model,
        profile,
        spine,
        &PipeLaw::RotationMinimizing,
        corners,
        1e-3,
        T,
    )
    .unwrap_or_else(|e| panic!("{corners:?}: {e}"))
    .shape
}

fn refusal(model: &mut Model, profile: &Shape, spine: &Shape, corners: PipeCorners) -> String {
    match make_pipe_shell_with(
        model,
        profile,
        spine,
        &PipeLaw::RotationMinimizing,
        corners,
        1e-3,
        T,
    ) {
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

fn close(got: f64, want: f64) -> bool {
    (got - want).abs() <= want.abs() * EXACT
}

/// An open polyline in the plane `z = 0`, square to a spine leaving the
/// origin along `+z`.
fn polyline(model: &mut Model, points: &[(f64, f64)]) -> Shape {
    let points: Vec<Point> = points.iter().map(|&(x, y)| Point::new(x, y, 0.0)).collect();
    make_polygon(model, &points, false, T).unwrap().shape
}

fn spine(model: &mut Model, points: &[Point]) -> Shape {
    make_polygon(model, points, false, T).unwrap().shape
}

/// A corner a spine turns: its angle, and where its outside lies in the
/// start section's own coordinates, the point's `α` being its distance on
/// that side of the axis.
#[derive(Clone, Copy)]
struct Corner {
    theta: f64,
    outward: (f64, f64),
}

impl Corner {
    fn alpha(self, (x, y): (f64, f64)) -> f64 {
        x * self.outward.0 + y * self.outward.1
    }

    /// The gain of one point's path.
    fn at(self, alpha: f64, round: bool) -> f64 {
        let k = 2.0 * (self.theta / 2.0).tan();
        if round {
            k * alpha.min(0.0) + self.theta * alpha.max(0.0)
        } else {
            k * alpha
        }
    }

    /// The gain integrated along a segment of length `s` over which `α`
    /// runs linearly from `a` to `b`.
    fn along(self, a: f64, b: f64, s: f64, round: bool) -> f64 {
        let k = 2.0 * (self.theta / 2.0).tan();
        let all = s * (a + b) / 2.0;
        if !round {
            return k * all;
        }
        let outer = if a >= 0.0 && b >= 0.0 {
            all
        } else if a <= 0.0 && b <= 0.0 {
            0.0
        } else {
            let (p, n) = if a > 0.0 { (a, b) } else { (b, a) };
            s * p * p / (2.0 * (p - n))
        };
        k * (all - outer) + self.theta * outer
    }
}

/// The closed forms of a polyline profile's sheet: its area and the
/// length of its boundary.
fn closed_form(profile: &[(f64, f64)], legs: f64, corners: &[Corner], round: bool) -> (f64, f64) {
    let mut length = 0.0;
    let mut area = 0.0;
    for pair in profile.windows(2) {
        let s = ((pair[1].0 - pair[0].0).powi(2) + (pair[1].1 - pair[0].1).powi(2)).sqrt();
        length += s;
        area += s * legs;
        for corner in corners {
            area += corner.along(corner.alpha(pair[0]), corner.alpha(pair[1]), s, round);
        }
    }
    let ends = [profile[0], profile[profile.len() - 1]];
    let rails: f64 = ends
        .iter()
        .map(|&e| legs + corners.iter().map(|c| c.at(c.alpha(e), round)).sum::<f64>())
        .sum();
    (area, 2.0 * length + rails)
}

const L: [(f64, f64); 3] = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0)];
const U: [(f64, f64); 4] = [(-1.0, 1.0), (-1.0, -1.0), (1.0, -1.0), (1.0, 1.0)];
const V: [(f64, f64); 3] = [(-1.0, 1.0), (0.0, -1.0), (1.0, 1.0)];
const SEGMENT: [(f64, f64); 2] = [(-1.0, 0.5), (1.5, 0.5)];

/// One corner, `l` up `+z` then `l` on, turned by `theta` toward `toward`
/// (a unit direction square to `z`): its outside lies away from `toward`.
fn one_corner(model: &mut Model, l: f64, theta: f64, toward: (f64, f64)) -> (Shape, Corner) {
    let p1 = Point::new(0.0, 0.0, l);
    let p2 = p1
        + Vector::new(
            l * theta.sin() * toward.0,
            l * theta.sin() * toward.1,
            l * theta.cos(),
        );
    let corner = Corner {
        theta,
        outward: (-toward.0, -toward.1),
    };
    (spine(model, &[Point::ORIGIN, p1, p2]), corner)
}

#[test]
fn an_open_profile_turns_one_corner_mitred_or_round_to_its_closed_form() {
    let l = 10.0;
    for profile in [&L[..], &U[..], &V[..], &SEGMENT[..]] {
        for theta in [FRAC_PI_2, FRAC_PI_3, 1e-3] {
            for (toward, reversed) in [
                ((1.0, 0.0), false),
                ((-1.0, 0.0), true),
                ((0.0, 1.0), false),
                ((0.6, -0.8), true),
            ] {
                for corners in [PipeCorners::Mitre, PipeCorners::Round] {
                    let mut model = Model::new();
                    let mut wire = polyline(&mut model, profile);
                    if reversed {
                        wire = wire.reversed();
                    }
                    let (path, corner) = one_corner(&mut model, l, theta, toward);
                    let sheet = sweep(&mut model, &wire, &path, corners);
                    let (area, boundary) = measured(&model, &sheet);
                    let round = corners == PipeCorners::Round;
                    let (want, around) = closed_form(profile, 2.0 * l, &[corner], round);
                    let case = format!("{profile:?} {theta} {toward:?} {reversed} {corners:?}");
                    assert!(close(area, want), "{case}: area {area} against {want}");
                    assert!(
                        close(boundary, around),
                        "{case}: boundary {boundary} against {around}"
                    );
                }
            }
        }
    }
}

/// A half circle of radius 1 from `(1, 0)` through `(0, 1)`, round a
/// right-angled corner toward `+x`: `α = -cos φ`, so mitred it gains
/// nothing (`∫ α = 0`) and rounded it loses `2` inside and gains `π/2`
/// outside.
#[test]
fn a_half_circle_turns_a_corner_to_its_closed_form() {
    for (corners, want) in [
        (PipeCorners::Mitre, 20.0 * PI),
        (PipeCorners::Round, 20.0 * PI - 2.0 + FRAC_PI_2),
    ] {
        let mut model = Model::new();
        let frame = ogeom_math::Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
        let circle = ogeom_math::Circle::new(frame, 1.0, T).unwrap();
        let edge = ogeom_algo::make_edge(
            &mut model,
            ogeom_geom::CircleCurve::new(circle).into(),
            (0.0, PI),
            T,
        )
        .unwrap()
        .shape;
        let wire = ogeom_algo::make_wire(&mut model, &[edge], T).unwrap().shape;
        let (path, _) = one_corner(&mut model, 10.0, FRAC_PI_2, (1.0, 0.0));
        let sheet = sweep(&mut model, &wire, &path, corners);
        let (area, _) = measured(&model, &sheet);
        assert!(close(area, want), "{corners:?}: {area} against {want}");
    }
}

/// Along one straight leg, or legs running straight on, there is no
/// corner: every way of turning one sweeps the prism of the profile, facing
/// as it does, a reversed profile facing the other way.
#[test]
fn an_open_profile_down_a_straight_spine_is_its_prism() {
    let mut cases: Vec<Vec<(f64, f64)>> = vec![L.to_vec(), SEGMENT.to_vec()];
    cases.push(vec![(-1.0, 0.0), (0.0, 0.0), (1.0, 0.0)]);
    for profile in &cases {
        for stops in [&[10.0][..], &[4.0, 10.0][..]] {
            for reversed in [false, true] {
                for corners in [
                    PipeCorners::Mitre,
                    PipeCorners::Round,
                    PipeCorners::Extended,
                ] {
                    let mut model = Model::new();
                    let mut wire = polyline(&mut model, profile);
                    if reversed {
                        wire = wire.reversed();
                    }
                    let mut points = vec![Point::ORIGIN];
                    points.extend(stops.iter().map(|&z| Point::new(0.0, 0.0, z)));
                    let path = spine(&mut model, &points);
                    let sheet = sweep(&mut model, &wire, &path, corners);
                    let prism = make_prism(&mut model, &wire, Vector::new(0.0, 0.0, 10.0), T)
                        .unwrap()
                        .shape;
                    let case = format!("{profile:?} {stops:?} {reversed} {corners:?}");
                    let (area, boundary) = measured(&model, &sheet);
                    let (want, around) = closed_form(profile, 10.0, &[], false);
                    assert!(close(area, want), "{case}: {area} against {want}");
                    assert!(
                        close(boundary, around),
                        "{case}: {boundary} against {around}"
                    );
                    // Each face of the sheet faces as the prism's face on
                    // the same plane does.
                    let facing = |shape: &Shape| -> Vec<(Point, Vector)> {
                        explore(&model, shape, Filter::OfType(ShapeType::Face))
                            .unwrap()
                            .iter()
                            .map(|f| ogeom_algo::face_normal(&model, f, T).unwrap())
                            .collect()
                    };
                    let reference = facing(&prism);
                    for (at, n) in facing(&sheet) {
                        let (_, m) = reference
                            .iter()
                            .find(|(p, m)| {
                                (at - *p).dot(*m).abs() < 1e-9 && n.cross(*m).magnitude() < 1e-9
                            })
                            .unwrap_or_else(|| panic!("{case}: no prism face on {at:?}"));
                        assert!(n.dot(*m) > 0.0, "{case}: a face turned over at {at:?}");
                    }
                }
            }
        }
    }
}

/// A profile running along the corner's axis (the line `x = 0` of the
/// section at the corner) and then out a flange inside the turn: the edge
/// on the axis turns on itself, and mitred or round the wall `x = 0` runs
/// its legs' `20` by `2`, the flange `y = 1` loses its triangle on each
/// leg: `40 + 2 · 9.5 = 59`.
#[test]
fn a_profile_edge_along_the_corner_axis_turns_on_itself() {
    let profile = [(0.0, -1.0), (0.0, 1.0), (1.0, 1.0)];
    for corners in [PipeCorners::Mitre, PipeCorners::Round] {
        let mut model = Model::new();
        let wire = polyline(&mut model, &profile);
        let (path, corner) = one_corner(&mut model, 10.0, FRAC_PI_2, (1.0, 0.0));
        let sheet = sweep(&mut model, &wire, &path, corners);
        let (area, boundary) = measured(&model, &sheet);
        let round = corners == PipeCorners::Round;
        let (want, around) = closed_form(&profile, 20.0, &[corner], round);
        assert!(close(want, 59.0), "{want}");
        assert!(close(area, want), "{corners:?}: {area} against {want}");
        assert!(
            close(boundary, around),
            "{corners:?}: {boundary} against {around}"
        );
    }
}

/// Two corners in one plane, turning the same way (a U) or back (a Z):
/// the section keeps its frame in the plane, so its outside at the second
/// corner is the first's for a U and the other side for a Z.
#[test]
fn an_open_profile_turns_two_corners_either_way() {
    let l = 10.0;
    for (name, far, second) in [
        ("U", Point::new(l, 0.0, 0.0), (-1.0, 0.0)),
        ("Z", Point::new(l, 0.0, 2.0 * l), (1.0, 0.0)),
    ] {
        for profile in [&L[..], &V[..]] {
            for corners in [PipeCorners::Mitre, PipeCorners::Round] {
                let mut model = Model::new();
                let wire = polyline(&mut model, profile);
                let path = spine(
                    &mut model,
                    &[
                        Point::ORIGIN,
                        Point::new(0.0, 0.0, l),
                        Point::new(l, 0.0, l),
                        far,
                    ],
                );
                let sheet = sweep(&mut model, &wire, &path, corners);
                let (area, boundary) = measured(&model, &sheet);
                let first = Corner {
                    theta: FRAC_PI_2,
                    outward: (-1.0, 0.0),
                };
                let then = Corner {
                    theta: FRAC_PI_2,
                    outward: second,
                };
                let round = corners == PipeCorners::Round;
                let (want, around) = closed_form(profile, 3.0 * l, &[first, then], round);
                let case = format!("{name} {profile:?} {corners:?}");
                assert!(close(area, want), "{case}: {area} against {want}");
                assert!(
                    close(boundary, around),
                    "{case}: {boundary} against {around}"
                );
            }
        }
    }
}

/// A 2 x 2 square about the spine cut into two open halves (two sides
/// each), along a spine turning corners out of any one plane.
///
/// Mitred, the square's points gain as much outside each corner as they
/// lose inside it, so the halves together are the perimeter times the
/// spine's length, `8 Σℓ`, whatever the frames do at the corners; the
/// square is turned against the corners to show it.
///
/// Rounded, along a staircase of right-angled corners about the
/// coordinate axes the square keeps its sides along the corners' axes: at
/// each corner the side outside turns a quarter drum `π`, the side inside
/// loses `2 · 2`, and each side across the turn loses its inner triangle
/// and gains its outer quarter disc, `π/4 - 1`: `8 Σℓ + 3π/2 - 6` a
/// corner. There every slice of the square across the turn reaches the
/// corner's axis, so the halves are also the walls of the closed square's
/// sweep, which the closed wire builds its own way, as the walls of its
/// face's solid.
#[test]
fn two_open_halves_make_the_closed_square_round_a_spine_out_of_plane() {
    let halves = |square: &[(f64, f64)]| {
        [
            vec![square[0], square[1], square[2]],
            vec![square[2], square[3], square[0]],
        ]
    };
    let length =
        |path: &[Point]| -> f64 { path.windows(2).map(|pair| pair[0].distance(pair[1])).sum() };
    let swept = |model: &mut Model, profile: &[(f64, f64)], line: &Shape, corners| {
        let wire = polyline(model, profile);
        let sheet = sweep(model, &wire, line, corners);
        measured(model, &sheet).0
    };
    {
        let turn = 0.4_f64;
        let square: Vec<(f64, f64)> = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
            .iter()
            .map(|&(x, y)| {
                (
                    x * turn.cos() - y * turn.sin(),
                    x * turn.sin() + y * turn.cos(),
                )
            })
            .collect();
        let path = [
            Point::ORIGIN,
            Point::new(0.0, 0.0, 10.0),
            Point::new(8.0, 0.0, 16.0),
            Point::new(8.0, 9.0, 16.0),
            Point::new(14.0, 12.0, 25.0),
        ];
        let mut model = Model::new();
        let line = spine(&mut model, &path);
        let area: f64 = halves(&square)
            .iter()
            .map(|half| swept(&mut model, half, &line, PipeCorners::Mitre))
            .sum();
        let want = 8.0 * length(&path);
        assert!(close(area, want), "mitred: {area} against {want}");
    }
    let square = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
    let path = [
        Point::ORIGIN,
        Point::new(0.0, 0.0, 10.0),
        Point::new(10.0, 0.0, 10.0),
        Point::new(10.0, 10.0, 10.0),
        Point::new(10.0, 10.0, 20.0),
    ];
    let mut model = Model::new();
    let line = spine(&mut model, &path);
    let area: f64 = halves(&square)
        .iter()
        .map(|half| swept(&mut model, half, &line, PipeCorners::Round))
        .sum();
    let want = 8.0 * length(&path) + 3.0 * (1.5 * PI - 6.0);
    assert!(close(area, want), "rounded: {area} against {want}");
    let closed = {
        let points: Vec<Point> = square.iter().map(|&(x, y)| Point::new(x, y, 0.0)).collect();
        make_polygon(&mut model, &points, true, T).unwrap().shape
    };
    let walls = sweep(&mut model, &closed, &line, PipeCorners::Round);
    assert_eq!(model.kind_of(&walls).unwrap(), ShapeType::Shell);
    let built = surface_properties(&model, &walls, Deflection::with_chord(1e-4).unwrap(), T)
        .unwrap()
        .mass;
    assert!(
        close(area, built),
        "rounded: {area} against the closed wire's {built}"
    );
}

/// The same sweep under a rigid motion, and down the same spine walked
/// from its other end (a wire reversed) where it is the motion of a
/// forward one: the closed form does not move.
#[test]
fn an_open_sweep_moved_or_walked_backward_keeps_its_closed_form() {
    let l = 10.0;
    let skew = Transform::translation(Vector::new(3.0, -7.0, 2.5))
        * Transform::rotation(
            Axis {
                location: Point::new(1.0, 2.0, 3.0),
                direction: Direction::new(Vector::new(1.0, 2.0, -0.5), T).unwrap(),
            },
            1.1,
        );
    // The spine (0,0,0) -> (0,0,l) -> (l,0,l) walked back starts at
    // (l,0,l) heading -x and turns to -z: the forward walk turned half
    // round the line through (l/2, 0, l/2) along (1, 0, -1), which takes
    // z to -x and x to -z and keeps the corner.
    let back = Transform::rotation(
        Axis {
            location: Point::new(l / 2.0, 0.0, l / 2.0),
            direction: Direction::new(Vector::new(1.0, 0.0, -1.0), T).unwrap(),
        },
        PI,
    );
    for corners in [PipeCorners::Mitre, PipeCorners::Round] {
        let round = corners == PipeCorners::Round;
        let corner = Corner {
            theta: FRAC_PI_2,
            outward: (-1.0, 0.0),
        };
        let (want, around) = closed_form(&L, 2.0 * l, &[corner], round);
        for (name, motion, walked_back) in [("moved", skew, false), ("walked back", back, true)] {
            let mut model = Model::new();
            let profile: Vec<Point> = L
                .iter()
                .map(|&(x, y)| motion.apply(Point::new(x, y, 0.0)))
                .collect();
            let wire = make_polygon(&mut model, &profile, false, T).unwrap().shape;
            let forward = [
                Point::ORIGIN,
                Point::new(0.0, 0.0, l),
                Point::new(l, 0.0, l),
            ];
            let path = if walked_back {
                spine(&mut model, &forward).reversed()
            } else {
                let moved: Vec<Point> = forward.iter().map(|p| motion.apply(*p)).collect();
                spine(&mut model, &moved)
            };
            let sheet = sweep(&mut model, &wire, &path, corners);
            let (area, boundary) = measured(&model, &sheet);
            let case = format!("{name} {corners:?}");
            assert!(close(area, want), "{case}: {area} against {want}");
            assert!(
                close(boundary, around),
                "{case}: {boundary} against {around}"
            );
        }
    }
}

/// A right-angled corner toward `+x` cuts the path of a point at `x = 1`
/// (inside the turn) short by 1 on each leg: a leg a hair longer still
/// sweeps, a leg a hair shorter is refused by name, either leg.
#[test]
fn a_leg_just_past_its_corners_reach_sweeps_and_one_short_of_it_is_refused() {
    let profile = [(-1.0, 0.0), (1.0, 0.0)];
    let corner = Corner {
        theta: FRAC_PI_2,
        outward: (-1.0, 0.0),
    };
    for corners in [PipeCorners::Mitre, PipeCorners::Round] {
        let round = corners == PipeCorners::Round;
        for (first, second) in [(1.001, 10.0), (10.0, 1.001), (0.999, 10.0), (10.0, 0.999)] {
            let mut model = Model::new();
            let wire = polyline(&mut model, &profile);
            let path = spine(
                &mut model,
                &[
                    Point::ORIGIN,
                    Point::new(0.0, 0.0, first),
                    Point::new(second, 0.0, first),
                ],
            );
            let case = format!("{corners:?} {first} {second}");
            if first.min(second) > 1.0 {
                let sheet = sweep(&mut model, &wire, &path, corners);
                let (area, _) = measured(&model, &sheet);
                let (want, _) = closed_form(&profile, first + second, &[corner], round);
                assert!(close(area, want), "{case}: {area} against {want}");
            } else {
                let message = refusal(&mut model, &wire, &path, corners);
                assert!(
                    message.contains("shorter than its corners' reach"),
                    "{case}: {message}"
                );
            }
        }
    }
}

/// A sweep refused by name: a profile's points, a spine's, the corners
/// asked for and the words of the refusal.
struct Refused<'a> {
    name: &'a str,
    profile: &'a [Point],
    path: &'a [Point],
    corners: PipeCorners,
    says: &'a str,
}

#[test]
fn open_wire_sweeps_it_cannot_build_are_refused_by_name() {
    let corner_path = [
        Point::ORIGIN,
        Point::new(0.0, 0.0, 10.0),
        Point::new(10.0, 0.0, 10.0),
    ];
    let across = [Point::new(-1.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0)];
    let cases = [
        Refused {
            name: "extended corner",
            profile: &across,
            path: &corner_path,
            corners: PipeCorners::Extended,
            says: "an extended corner squares the outside of a solid's section; an open \
                   wire sweeps a sheet, whose corners are mitred or round",
        },
        Refused {
            name: "leaning segment",
            profile: &[Point::new(-1.0, 0.0, -0.1), Point::new(1.0, 0.0, 0.1)],
            path: &corner_path,
            corners: PipeCorners::Mitre,
            says: "the profile leans along its spine",
        },
        Refused {
            name: "segment off the start",
            profile: &[Point::new(-1.0, 0.0, 1.0), Point::new(1.0, 0.0, 1.0)],
            path: &corner_path,
            corners: PipeCorners::Round,
            says: "the profile does not sit at the spine's start",
        },
        Refused {
            name: "doubled back",
            profile: &across,
            path: &[
                Point::ORIGIN,
                Point::new(0.0, 0.0, 10.0),
                Point::new(1.0, 0.0, 0.0),
            ],
            corners: PipeCorners::Mitre,
            says: "the spine all but doubles back at a corner",
        },
        // The last leg runs back at z = 5 through the first leg's wall on
        // x = 0.
        Refused {
            name: "winding back through itself",
            profile: &[Point::new(0.0, -2.0, 0.0), Point::new(0.0, 2.0, 0.0)],
            path: &[
                Point::ORIGIN,
                Point::new(0.0, 0.0, 10.0),
                Point::new(10.0, 0.0, 10.0),
                Point::new(10.0, 0.0, 5.0),
                Point::new(-5.0, 0.0, 5.0),
            ],
            corners: PipeCorners::Mitre,
            says: "the swept sheet runs into itself",
        },
    ];
    for case in cases {
        let mut model = Model::new();
        let wire = make_polygon(&mut model, case.profile, false, T)
            .unwrap()
            .shape;
        let line = spine(&mut model, case.path);
        let message = refusal(&mut model, &wire, &line, case.corners);
        assert!(message.contains(case.says), "{}: {message}", case.name);
    }
}

/// A closed wire walked backward is still closed: mitred down an L it
/// sweeps what it sweeps walked forward, not an open wire's sheet.
#[test]
fn a_closed_wire_walked_backward_is_not_taken_for_an_open_one() {
    let square = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
    let mut areas = Vec::new();
    for reversed in [false, true] {
        let mut model = Model::new();
        let points: Vec<Point> = square.iter().map(|&(x, y)| Point::new(x, y, 0.0)).collect();
        let mut wire = make_polygon(&mut model, &points, true, T).unwrap().shape;
        if reversed {
            wire = wire.reversed();
        }
        let (path, _) = one_corner(&mut model, 10.0, FRAC_PI_2, (1.0, 0.0));
        let swept = sweep(&mut model, &wire, &path, PipeCorners::Mitre);
        areas.push(
            surface_properties(&model, &swept, Deflection::with_chord(1e-4).unwrap(), T)
                .unwrap()
                .mass,
        );
    }
    assert!(close(areas[1], areas[0]), "{areas:?}");
}

/// A segment lying in the plane the spine turns in sweeps a flat sheet,
/// rounds and all: every face of it faces as the segment's prism down the
/// first leg does, and the other way for the segment reversed.
#[test]
fn an_open_sheet_faces_as_its_prism_round_its_corners() {
    let l = 10.0;
    for reversed in [false, true] {
        for corners in [PipeCorners::Mitre, PipeCorners::Round] {
            let mut model = Model::new();
            let mut wire = polyline(&mut model, &[(-1.0, 0.0), (1.0, 0.0)]);
            if reversed {
                wire = wire.reversed();
            }
            let path = spine(
                &mut model,
                &[
                    Point::ORIGIN,
                    Point::new(0.0, 0.0, l),
                    Point::new(l, 0.0, l),
                    Point::new(l, 0.0, 0.0),
                ],
            );
            let sheet = sweep(&mut model, &wire, &path, corners);
            let prism = make_prism(&mut model, &wire, Vector::new(0.0, 0.0, l), T)
                .unwrap()
                .shape;
            let face = explore(&model, &prism, Filter::OfType(ShapeType::Face)).unwrap();
            let (_, want) = ogeom_algo::face_normal(&model, &face[0], T).unwrap();
            let faces = explore(&model, &sheet, Filter::OfType(ShapeType::Face)).unwrap();
            assert!(faces.len() >= 3, "{reversed} {corners:?}");
            for face in &faces {
                let (at, n) = ogeom_algo::face_normal(&model, face, T).unwrap();
                assert!(
                    n.cross(want).magnitude() < 1e-9 && n.dot(want) > 0.0,
                    "{reversed} {corners:?}: {n:?} at {at:?} against {want:?}"
                );
            }
        }
    }
}

/// Every vertex of a shape where it stands.
fn vertex_points(model: &Model, shape: &Shape) -> Vec<Point> {
    explore(model, shape, Filter::OfType(ShapeType::Vertex))
        .unwrap()
        .iter()
        .map(|v| {
            let at = model.node(v).unwrap().data().as_vertex().unwrap().point;
            v.transform(model.datums()).unwrap().apply(at)
        })
        .collect()
}

/// A spine or profile placed by a location, its own or its edges' under
/// its own, sweeps what the same shapes built where they stand sweep: the
/// same area, boundary and vertices.
#[test]
fn an_open_sweep_of_located_shapes_is_the_sweep_where_they_stand() {
    let l = 10.0;
    let forward = [
        Point::ORIGIN,
        Point::new(0.0, 0.0, l),
        Point::new(l, 0.0, l),
    ];
    let shift = Transform::translation(Vector::new(5.0, 0.0, 0.0));
    let turn = Transform::rotation(
        Axis {
            location: Point::new(1.0, 2.0, 3.0),
            direction: Direction::new(Vector::new(1.0, 2.0, -0.5), T).unwrap(),
        },
        1.1,
    );
    // Which of the profile and the spine carries the placement, and how:
    // on the wire, on the wire over one on each edge, or on the wire
    // over another on the wire.
    #[derive(Clone, Copy, Debug)]
    enum Placed {
        Not,
        Wire,
        WireOverEdges,
        WireOverWire,
    }
    // Each placement is an outer motion over an inner one.
    let cases = [
        (
            "spine shifted",
            shift,
            Transform::IDENTITY,
            Placed::Not,
            Placed::Wire,
        ),
        (
            "spine turned",
            turn,
            Transform::IDENTITY,
            Placed::Not,
            Placed::Wire,
        ),
        (
            "spine over its edges",
            turn,
            shift,
            Placed::Not,
            Placed::WireOverEdges,
        ),
        (
            "profile turned",
            turn,
            Transform::IDENTITY,
            Placed::Wire,
            Placed::Not,
        ),
        (
            "both, nested",
            shift,
            turn,
            Placed::WireOverWire,
            Placed::WireOverEdges,
        ),
    ];
    // A shape of `points` placed by `outer` over `inner` as `how` asks:
    // baked in where `how` places nothing.
    let build =
        |model: &mut Model, points: &[Point], outer: Transform, inner: Transform, how: Placed| {
            let (outer, inner) = match how {
                Placed::Not => {
                    let moved: Vec<Point> =
                        points.iter().map(|p| (outer * inner).apply(*p)).collect();
                    return make_polygon(model, &moved, false, T).unwrap().shape;
                }
                Placed::Wire => (outer * inner, Transform::IDENTITY),
                Placed::WireOverEdges | Placed::WireOverWire => (outer, inner),
            };
            let (outer, inner) = (
                Location::of(model.add_datum(outer)),
                Location::of(model.add_datum(inner)),
            );
            let wire = make_polygon(model, points, false, T).unwrap().shape;
            let wire = match how {
                Placed::WireOverEdges => {
                    let edges: Vec<Shape> = model
                        .ordered_children_of(&wire)
                        .unwrap()
                        .iter()
                        .map(|e| e.moved(&inner))
                        .collect();
                    make_wire(model, &edges, T).unwrap().shape
                }
                Placed::WireOverWire => wire.moved(&inner),
                _ => wire,
            };
            wire.moved(&outer)
        };
    for corners in [PipeCorners::Mitre, PipeCorners::Round] {
        let round = corners == PipeCorners::Round;
        let corner = Corner {
            theta: FRAC_PI_2,
            outward: (-1.0, 0.0),
        };
        let (want, around) = closed_form(&L, 2.0 * l, &[corner], round);
        let profile: Vec<Point> = L.iter().map(|&(x, y)| Point::new(x, y, 0.0)).collect();
        for (name, outer, inner, profile_placed, spine_placed) in cases {
            let case = format!("{name} {corners:?}");
            let mut model = Model::new();
            let wire = build(&mut model, &profile, outer, inner, profile_placed);
            let path = build(&mut model, &forward, outer, inner, spine_placed);
            let sheet = sweep(&mut model, &wire, &path, corners);
            let (area, boundary) = measured(&model, &sheet);
            assert!(close(area, want), "{case}: {area} against {want}");
            assert!(
                close(boundary, around),
                "{case}: {boundary} against {around}"
            );
            let mut baked = Model::new();
            let wire = build(&mut baked, &profile, outer, inner, Placed::Not);
            let path = build(&mut baked, &forward, outer, inner, Placed::Not);
            let where_they_stand = sweep(&mut baked, &wire, &path, corners);
            let got = vertex_points(&model, &sheet);
            let expected = vertex_points(&baked, &where_they_stand);
            assert_eq!(got.len(), expected.len(), "{case}");
            for p in &expected {
                assert!(
                    got.iter().any(|q| q.distance(*p) <= 1e-9),
                    "{case}: no vertex at {p:?}"
                );
            }
        }
    }
}

/// The L inside the turn, touching the corner's axis at one end.
const INSIDE: [(f64, f64); 3] = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)];

/// A corner turning by next to nothing, as far down as an angle's
/// tolerance: mitred it sweeps to its closed form, every point's path held
/// to the confusion distance, so the area to the profile's length times
/// twice that and the boundary to four times it; rounded, where the round's
/// arc or its cut inside the turn falls within the confusion distance of
/// the square end, it is refused by name. A turn below an angle's
/// tolerance is none, and a turn of `1e-4` rounds to its closed form.
#[test]
fn a_corner_turning_next_to_nothing_mitres_or_is_refused_by_name() {
    let l = 10.0;
    for profile in [&L[..], &SEGMENT[..], &INSIDE[..]] {
        let length: f64 = profile
            .windows(2)
            .map(|pair| ((pair[1].0 - pair[0].0).powi(2) + (pair[1].1 - pair[0].1).powi(2)).sqrt())
            .sum();
        for theta in [1e-13, 1e-11, 1e-9, 1e-8, 1e-7, 1e-4] {
            for corners in [PipeCorners::Mitre, PipeCorners::Round] {
                let case = format!("{profile:?} {theta} {corners:?}");
                let mut model = Model::new();
                let wire = polyline(&mut model, profile);
                let (path, corner) = one_corner(&mut model, l, theta, (1.0, 0.0));
                let round = corners == PipeCorners::Round;
                if round && (1e-11..=1e-8).contains(&theta) {
                    let message = refusal(&mut model, &wire, &path, corners);
                    assert!(
                        message.contains("a round corner turns so little"),
                        "{case}: {message}"
                    );
                    continue;
                }
                // A round of 1e-7 at a distance of 1 sits on the confusion
                // distance itself.
                if round && theta == 1e-7 {
                    continue;
                }
                let sheet = sweep(&mut model, &wire, &path, corners);
                let (area, boundary) = measured(&model, &sheet);
                let (want, around) = closed_form(profile, 2.0 * l, &[corner], round);
                let exact = theta == 1e-13 || theta == 1e-4;
                let (area_slack, boundary_slack) = if exact {
                    (want * EXACT, around * EXACT)
                } else {
                    (length * 2.0 * T.confusion(), 4.0 * T.confusion())
                };
                assert!(
                    (area - want).abs() <= area_slack,
                    "{case}: area {area} against {want}"
                );
                assert!(
                    (boundary - around).abs() <= boundary_slack,
                    "{case}: boundary {boundary} against {around}"
                );
            }
        }
    }
}

/// A piece of a spine: a line to a point, or an arc of a radius through an
/// angle, bending toward a direction, after turning sharply to a heading
/// where one is given.
#[derive(Clone, Copy)]
enum Piece {
    Line(Point),
    Arc {
        turned: Option<Vector>,
        toward: Vector,
        radius: f64,
        angle: f64,
    },
}

/// A spine from the origin heading `+z`, piece by piece, each starting
/// where the last ends.
fn pieced(model: &mut Model, pieces: &[Piece]) -> Shape {
    let mut at = Point::ORIGIN;
    let mut heading = Vector::new(0.0, 0.0, 1.0);
    let mut from = ogeom_algo::make_vertex(model, at).shape;
    let mut edges = Vec::new();
    for piece in pieces {
        let (curve, range, to): (ogeom_geom::Curve, _, _) = match *piece {
            Piece::Line(to) => {
                let line = ogeom_geom::LineCurve::segment(at, to, T).unwrap();
                heading = (to - at) / to.distance(at);
                let range = line.domain();
                (line.into(), range, to)
            }
            Piece::Arc {
                turned,
                toward,
                radius,
                angle,
            } => {
                if let Some(h) = turned {
                    heading = h;
                }
                let centre = at + toward * radius;
                let x = (at - centre) / radius;
                let frame = ogeom_math::Frame::new(
                    centre,
                    Direction::new(x.cross(heading), T).unwrap(),
                    Direction::new(x, T).unwrap(),
                    T,
                )
                .unwrap();
                let circle = ogeom_math::Circle::new(frame, radius, T).unwrap();
                let curve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(circle).into();
                let to = curve.point_at(angle, T).unwrap();
                let d = curve.d1_at(angle, T).unwrap();
                heading = d / d.magnitude();
                (curve, (0.0, angle), to)
            }
        };
        let end = ogeom_algo::make_vertex(model, to).shape;
        edges.push(
            ogeom_algo::make_edge_between(model, curve, range, &from, &end, T)
                .unwrap()
                .shape,
        );
        at = to;
        from = end;
    }
    make_wire(model, &edges, T).unwrap().shape
}

const X: Vector = Vector::new(1.0, 0.0, 0.0);
const Y: Vector = Vector::new(0.0, 1.0, 0.0);
const Z: Vector = Vector::new(0.0, 0.0, 1.0);

/// A profile wholly outside a right-angled corner toward `+x` (`α = -x`),
/// and its mirror, wholly outside one toward `-x` (`α = x`).
const OUTSIDE: [(f64, f64); 3] = [(-3.0, -1.0), (-1.0, -1.0), (-1.0, 1.0)];
const MIRRORED: [(f64, f64); 3] = [(3.0, -1.0), (1.0, -1.0), (1.0, 1.0)];

/// The area and boundary of a polyline profile's sheet where every point
/// sweeps a path whose length is linear in the point, `run(x, y)`: each
/// segment sweeps its length times the run of its middle, and the
/// boundary is the profile twice and the runs of its two ends.
fn linear_closed_form(profile: &[(f64, f64)], run: impl Fn(f64, f64) -> f64) -> (f64, f64) {
    let mut length = 0.0;
    let mut area = 0.0;
    for pair in profile.windows(2) {
        let s = ((pair[1].0 - pair[0].0).powi(2) + (pair[1].1 - pair[0].1).powi(2)).sqrt();
        length += s;
        area += s * run(
            f64::midpoint(pair[0].0, pair[1].0),
            f64::midpoint(pair[0].1, pair[1].1),
        );
    }
    let ends = [profile[0], profile[profile.len() - 1]];
    (
        area,
        2.0 * length + run(ends[0].0, ends[0].1) + run(ends[1].0, ends[1].1),
    )
}

/// Down legs of lines and arcs, a point's path is its run down each line,
/// its arc about each arc's axis (the angle times its distance from the
/// axis, linear in the point while the section stays on one side of it)
/// and, at each right-angled corner it lies outside, `2 α` mitred
/// (`α tan(θ/2)` on each side) or `π α / 2` rounded.
///
/// Up `+z` 10, then sharply `+x` into an arc of radius 8 bending toward
/// `+z`, `-z` or `+y`: the corner turns the section's `(x, y)` to the
/// offset `(0, y, -x)` from it, `α = -x`, at `8 + x`, `8 - x` or `8 - y`
/// from the arc's axis. An arc of radius 8 toward `+x` first takes it to
/// `(0, y, -x)` at `8 - x` from the axis; a sharp turn up `+z` then has
/// `α = x` and takes it to `(x, y, 0)`, and an arc of radius 5 bending
/// toward `-x` from there holds it at `5 + x`. Up `+z`, sharply `+x` into a
/// quarter arc toward `+y` (at `8 - y`) and sharply down `-z` turns two
/// corners with `α = -x`. Tangent, a line and an arc turn no corner.
#[test]
fn an_open_profile_outside_corners_against_arcs_sweeps_to_its_closed_form() {
    let l = 10.0;
    let q = FRAC_PI_2;
    let arc = |turned, toward, radius, angle| Piece::Arc {
        turned,
        toward,
        radius,
        angle,
    };
    let up = Point::new(0.0, 0.0, l);
    let gain = |alpha: f64, round: bool| if round { q * alpha } else { 2.0 * alpha };
    // A case's name, profile, spine, and the run of a point `(x, y)`
    // gaining `g` round the corner.
    type Case<'a> = (&'a str, &'a [(f64, f64)], Vec<Piece>, Run);
    type Run = fn(f64, f64, f64) -> f64;
    let cases: Vec<Case<'_>> = vec![
        (
            "line, arc up",
            &OUTSIDE,
            vec![Piece::Line(up), arc(Some(X), Z, 8.0, q)],
            |x, _, g| 10.0 + g + FRAC_PI_2 * (8.0 + x),
        ),
        (
            "line, arc down",
            &OUTSIDE,
            vec![Piece::Line(up), arc(Some(X), -Z, 8.0, 0.5)],
            |x, _, g| 10.0 + g + 0.5 * (8.0 - x),
        ),
        (
            "line, arc aside",
            &OUTSIDE,
            vec![Piece::Line(up), arc(Some(X), Y, 8.0, q)],
            |_, y, g| 10.0 + g + FRAC_PI_2 * (8.0 - y),
        ),
        (
            "arc, line",
            &MIRRORED,
            vec![
                arc(None, X, 8.0, q),
                Piece::Line(Point::new(8.0, 0.0, 8.0 + l)),
            ],
            |x, _, g| FRAC_PI_2 * (8.0 - x) + g + 10.0,
        ),
        (
            "arc, arc",
            &MIRRORED,
            vec![arc(None, X, 8.0, q), arc(Some(Z), -X, 5.0, q)],
            |x, _, g| FRAC_PI_2 * (8.0 - x) + g + FRAC_PI_2 * (5.0 + x),
        ),
    ];
    for (name, profile, pieces, run) in &cases {
        for corners in [PipeCorners::Mitre, PipeCorners::Round] {
            let round = corners == PipeCorners::Round;
            let mut model = Model::new();
            let wire = polyline(&mut model, profile);
            let path = pieced(&mut model, pieces);
            let sheet = sweep(&mut model, &wire, &path, corners);
            let (area, boundary) = measured(&model, &sheet);
            let alpha_sign = if profile[0].0 < 0.0 { -1.0 } else { 1.0 };
            let (want, around) =
                linear_closed_form(profile, |x, y| run(x, y, gain(alpha_sign * x, round)));
            let case = format!("{name} {corners:?}");
            assert!(close(area, want), "{case}: area {area} against {want}");
            assert!(
                close(boundary, around),
                "{case}: boundary {boundary} against {around}"
            );
        }
    }
    // Two corners, both with `α = -x`.
    for corners in [PipeCorners::Mitre, PipeCorners::Round] {
        let round = corners == PipeCorners::Round;
        let mut model = Model::new();
        let wire = polyline(&mut model, &OUTSIDE);
        let path = pieced(
            &mut model,
            &[
                Piece::Line(up),
                arc(Some(X), Y, 8.0, q),
                Piece::Line(Point::new(8.0, 8.0, 0.0)),
            ],
        );
        // The second corner starts the last line straight down.
        let sheet = sweep(&mut model, &wire, &path, corners);
        let (area, boundary) = measured(&model, &sheet);
        let (want, around) = linear_closed_form(&OUTSIDE, |x, y| {
            2.0 * l + 2.0 * gain(-x, round) + q * (8.0 - y)
        });
        assert!(close(area, want), "{corners:?}: area {area} against {want}");
        assert!(
            close(boundary, around),
            "{corners:?}: boundary {boundary} against {around}"
        );
    }
    // Tangent: a line into an arc toward `+x` (at `8 - x`) and an arc into
    // a line, the profile on both sides of the spine.
    for (pieces, run) in [
        (
            vec![Piece::Line(up), arc(None, X, 8.0, q)],
            (|x, _| 10.0 + FRAC_PI_2 * (8.0 - x)) as fn(f64, f64) -> f64,
        ),
        (
            vec![
                arc(None, X, 8.0, q),
                Piece::Line(Point::new(8.0 + l, 0.0, 8.0)),
            ],
            |x, _| FRAC_PI_2 * (8.0 - x) + 10.0,
        ),
    ] {
        for corners in [PipeCorners::Mitre, PipeCorners::Round] {
            let mut model = Model::new();
            let wire = polyline(&mut model, &L);
            let path = pieced(&mut model, &pieces);
            let sheet = sweep(&mut model, &wire, &path, corners);
            let (area, boundary) = measured(&model, &sheet);
            let (want, around) = linear_closed_form(&L, run);
            assert!(close(area, want), "{corners:?}: area {area} against {want}");
            assert!(
                close(boundary, around),
                "{corners:?}: boundary {boundary} against {around}"
            );
        }
    }
}

/// A rectangle wholly outside a corner against an arc, swept as a face,
/// is bounded by its two ends and the walls its two open halves sweep:
/// rounded, exactly; mitred, within what the face's skin is held to.
#[test]
fn two_open_halves_outside_a_corner_against_an_arc_bound_the_face_swept() {
    let rectangle = [(-3.0, -1.0), (-1.0, -1.0), (-1.0, 1.0), (-3.0, 1.0)];
    let halves = [
        vec![rectangle[0], rectangle[1], rectangle[2]],
        vec![rectangle[2], rectangle[3], rectangle[0]],
    ];
    let pieces = |toward: Vector| {
        [
            Piece::Line(Point::new(0.0, 0.0, 10.0)),
            Piece::Arc {
                turned: Some(X),
                toward,
                radius: 8.0,
                angle: FRAC_PI_2,
            },
        ]
    };
    for toward in [Z, Y] {
        for corners in [PipeCorners::Mitre, PipeCorners::Round] {
            let mut model = Model::new();
            let path = pieced(&mut model, &pieces(toward));
            let sheets: f64 = halves
                .iter()
                .map(|half| {
                    let wire = polyline(&mut model, half);
                    let sheet = sweep(&mut model, &wire, &path, corners);
                    measured(&model, &sheet).0
                })
                .sum();
            let face = {
                let points: Vec<Point> = rectangle
                    .iter()
                    .map(|&(x, y)| Point::new(x, y, 0.0))
                    .collect();
                let wire = make_polygon(&mut model, &points, true, T).unwrap().shape;
                let plane = ogeom_math::Plane::new(
                    ogeom_math::Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap(),
                );
                ogeom_algo::make_face(
                    &mut model,
                    ogeom_geom::PlaneSurface::new(plane).into(),
                    &[wire],
                    T,
                )
                .unwrap()
                .shape
            };
            let solid = sweep(&mut model, &face, &path, corners);
            let walls =
                surface_properties(&model, &solid, Deflection::with_chord(1e-4).unwrap(), T)
                    .unwrap()
                    .mass
                    - 2.0 * 4.0;
            // The skin stands within its tolerance of the walls, whose
            // boundary is the rectangle's perimeter twice and the paths of
            // its four corners, each under 30 long.
            let bound = if corners == PipeCorners::Mitre {
                1e-3 * (2.0 * 8.0 + 4.0 * 30.0)
            } else {
                walls * EXACT
            };
            assert!(
                (sheets - walls).abs() <= bound,
                "{toward:?} {corners:?}: halves {sheets} against the face's walls {walls}"
            );
        }
    }
}

/// `f` integrated over `[a, b]` by 5-point Gauss–Legendre on 256 panels:
/// for the analytic path lengths below, exact to rounding.
fn integral(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    const NODES: [(f64, f64); 5] = [
        (0.0, 0.568_888_888_888_888_9),
        (0.538_469_310_105_683, 0.478_628_670_499_366_5),
        (-0.538_469_310_105_683, 0.478_628_670_499_366_5),
        (0.906_179_845_938_664, 0.236_926_885_056_189_1),
        (-0.906_179_845_938_664, 0.236_926_885_056_189_1),
    ];
    let panels = 256;
    let h = (b - a) / f64::from(panels);
    let mut sum = 0.0;
    for k in 0..panels {
        let mid = a + h * (f64::from(k) + 0.5);
        for (x, w) in NODES {
            sum += w * f(mid + x * h / 2.0);
        }
    }
    sum * h / 2.0
}

/// The area and boundary of a polyline profile's sheet whose points sweep
/// paths of length `run(x, y)`: each segment, split where it crosses
/// `x = 0` (where the path's length has a kink), sweeps the integral of
/// the run along it, and the boundary is the profile twice and the runs
/// of its two ends.
fn integrated_closed_form(profile: &[(f64, f64)], run: impl Fn(f64, f64) -> f64) -> (f64, f64) {
    let mut length = 0.0;
    let mut area = 0.0;
    for pair in profile.windows(2) {
        let ((x0, y0), (x1, y1)) = (pair[0], pair[1]);
        let s = (x1 - x0).hypot(y1 - y0);
        length += s;
        let at = |t: f64| run(x0 + (x1 - x0) * t, y0 + (y1 - y0) * t);
        let mut breaks = vec![0.0, 1.0];
        if x0 * x1 < 0.0 {
            breaks.insert(1, x0 / (x0 - x1));
        }
        for w in breaks.windows(2) {
            area += s * integral(at, w[0], w[1]);
        }
    }
    let ends = [profile[0], profile[profile.len() - 1]];
    (
        area,
        2.0 * length + run(ends[0].0, ends[0].1) + run(ends[1].0, ends[1].1),
    )
}

/// Inside a corner against an arc turning in the arc's plane, each leg's
/// wall runs on until it meets the other's: a point's two paths end where
/// they cross. The spines are those of the closed-form test outside such
/// corners, and a point inside the corner at `s` from the axis crosses at
/// these closed forms:
///
/// - up `+z` 10 and sharply `+x` into an arc of radius 8 toward `+z`: the
///   line `x = s` meets the circle of radius `8 + s` about `(0, 0, 18)` at
///   `z = 18 - √(64 + 16s)`, the arc there turned `atan(s / √(64 + 16s))`;
/// - the same into an arc toward `-z` of `0.5` rad: the circle of radius
///   `8 - s` about `(0, 0, 2)`, at `z = 2 + √(64 - 16s)`, turned
///   `atan(s / √(64 - 16s))`;
/// - a quarter arc of radius 8 toward `+x`, then sharply up `+z` 10: the
///   circle of radius `8 + s` about `(8, 0, 0)` meets the line `x = 8 - s`
///   at `z = √(64 + 16s)`, turned `atan2(√(64 + 16s), s)`;
/// - the same into a quarter arc of radius 5 toward `-x`: the circles of
///   radius `8 + s` about `(8, 0, 0)` and `5 - s` about `(3, 0, 8)`.
///
/// The profiles reach both sides of the corner's axis, along it and across
/// it (a parabola or an ellipse, in or across the corner's plane), and
/// points outside run on as the closed-form test's: mitred or round the
/// sheet is that integral.
/// A point's path round a corner against an arc, by its `x` and `y` in the
/// start section and whether the corner is rounded.
type Run = Box<dyn Fn(f64, f64, bool) -> f64>;

/// The spines of the closed-form test inside corners against arcs, each
/// with whether its corner turns toward `-x` (the profile mirrored in `x`
/// to reach inside it) and the closed-form run of a point.
fn inside_cases() -> Vec<(&'static str, Vec<Piece>, bool, Run)> {
    let q = FRAC_PI_2;
    let arc = |turned, toward, radius, angle| Piece::Arc {
        turned,
        toward,
        radius,
        angle,
    };
    let up = Point::new(0.0, 0.0, 10.0);
    // The run of a point outside the corner at `α`: mitred `2α`, round
    // `π α / 2`.
    let gain = move |alpha: f64, round: bool| if round { q * alpha } else { 2.0 * alpha };
    // The two circles' crossing nearest the corner `(8, 8)` in `(x, z)`.
    let circles = |s: f64| {
        let (c1, r1, c2, r2) = ((8.0, 0.0), 8.0 + s, (3.0, 8.0), 5.0 - s);
        let (ex, ez) = (c2.0 - c1.0, c2.1 - c1.1);
        let d = f64::hypot(ex, ez);
        let a = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
        let h = (r1 * r1 - a * a).sqrt();
        let (ux, uz) = (ex / d, ez / d);
        // Of the two, the one on the corner's side of the centres' line.
        (c1.0 + a * ux + h * uz, c1.1 + a * uz - h * ux)
    };
    vec![
        (
            "line, arc up",
            vec![Piece::Line(up), arc(Some(X), Z, 8.0, q)],
            false,
            Box::new(move |x, _, round| {
                if x <= 0.0 {
                    return 10.0 + gain(-x, round) + q * (8.0 + x);
                }
                let root = (64.0 + 16.0 * x).sqrt();
                18.0 - root + (8.0 + x) * (q - (x / root).atan())
            }),
        ),
        (
            "line, arc down",
            vec![Piece::Line(up), arc(Some(X), -Z, 8.0, 0.5)],
            false,
            Box::new(move |x, _, round| {
                if x <= 0.0 {
                    return 10.0 + gain(-x, round) + 0.5 * (8.0 - x);
                }
                let root = (64.0 - 16.0 * x).sqrt();
                2.0 + root + (8.0 - x) * (0.5 - (x / root).atan())
            }),
        ),
        (
            "arc, line",
            vec![
                arc(None, X, 8.0, q),
                Piece::Line(Point::new(8.0, 0.0, 18.0)),
            ],
            true,
            Box::new(move |x, _, round| {
                if x >= 0.0 {
                    return q * (8.0 - x) + gain(x, round) + 10.0;
                }
                let (s, root) = (-x, (64.0 - 16.0 * x).sqrt());
                (8.0 + s) * root.atan2(s) + 18.0 - root
            }),
        ),
        (
            "arc, arc",
            vec![arc(None, X, 8.0, q), arc(Some(Z), -X, 5.0, q)],
            true,
            Box::new(move |x, _, round| {
                if x >= 0.0 {
                    return q * (8.0 - x) + gain(x, round) + q * (5.0 + x);
                }
                let s = -x;
                let (cx, cz) = circles(s);
                (8.0 + s) * cz.atan2(8.0 - cx) + (5.0 - s) * (q - (cz - 8.0).atan2(cx - 3.0))
            }),
        ),
        (
            "line, arc up, line on",
            vec![
                Piece::Line(up),
                arc(Some(X), Z, 8.0, q),
                Piece::Line(Point::new(18.0, 0.0, 18.0)),
            ],
            false,
            Box::new(move |x, _, round| {
                if x <= 0.0 {
                    return 20.0 + 2.0 * gain(-x, round) + q * (8.0 + x);
                }
                let root = (64.0 + 16.0 * x).sqrt();
                2.0 * (18.0 - root) + (8.0 + x) * (q - 2.0 * (x / root).atan())
            }),
        ),
        (
            "line, arc up, line back",
            vec![
                Piece::Line(up),
                arc(Some(X), Z, 8.0, q),
                Piece::Line(Point::new(-2.0, 0.0, 18.0)),
            ],
            false,
            Box::new(move |x, _, round| {
                if x >= 0.0 {
                    let root = (64.0 + 16.0 * x).sqrt();
                    return 18.0 - root
                        + (8.0 + x) * (q - (x / root).atan())
                        + gain(x, round)
                        + 10.0;
                }
                let (s, root) = (-x, (64.0 + 16.0 * x).sqrt());
                10.0 + gain(s, round) + (8.0 - s) * root.atan2(s) + root + 2.0
            }),
        ),
    ]
}

#[test]
fn an_open_profile_inside_a_corner_against_an_arc_sweeps_to_its_closed_form() {
    let inside = [(-1.0, -1.0), (2.0, -1.0), (2.0, 1.0)];
    let across = [(-1.0, 1.5), (2.5, -0.5)];
    let mirrored = |p: &[(f64, f64)]| p.iter().map(|&(x, y)| (-x, y)).collect::<Vec<_>>();
    let cases = inside_cases();
    for (name, pieces, flipped, run) in &cases {
        for profile in [&inside[..], &across[..]] {
            let profile = if *flipped {
                mirrored(profile)
            } else {
                profile.to_vec()
            };
            for corners in [PipeCorners::Mitre, PipeCorners::Round] {
                let round = corners == PipeCorners::Round;
                let mut model = Model::new();
                let wire = polyline(&mut model, &profile);
                let path = pieced(&mut model, pieces);
                let sheet = sweep(&mut model, &wire, &path, corners);
                let (area, boundary) = measured(&model, &sheet);
                let (want, around) = integrated_closed_form(&profile, |x, y| run(x, y, round));
                let case = format!("{name} {profile:?} {corners:?}");
                assert!(close(area, want), "{case}: area {area} against {want}");
                assert!(
                    close(boundary, around),
                    "{case}: boundary {boundary} against {around}"
                );
            }
        }
    }
}

/// A curved edge of the profile in the plane `z = 0`: an arc of a circle
/// about `(x, y)` of a radius from one angle to another, or a cubic Bézier
/// through four control points.
#[derive(Clone, Copy)]
enum Bend {
    Arc((f64, f64), f64, (f64, f64)),
    Bezier([(f64, f64); 4]),
}

impl Bend {
    /// The bend mirrored in `x`.
    fn mirrored(self) -> Self {
        match self {
            Bend::Arc((x, y), r, angles) => Bend::Arc((-x, y), r, angles),
            Bend::Bezier(points) => Bend::Bezier(points.map(|(x, y)| (-x, y))),
        }
    }

    /// The point at `t` in `[0, 1]` and the speed there, an arc's `x` about
    /// its centre taken with the sign `mirror`.
    fn at(self, t: f64, mirror: f64) -> ((f64, f64), f64) {
        match self {
            Bend::Arc((x, y), r, (a, b)) => {
                let phi = a + (b - a) * t;
                (
                    (x + mirror * r * phi.cos(), y + r * phi.sin()),
                    r * (b - a).abs(),
                )
            }
            Bend::Bezier(p) => {
                let u = 1.0 - t;
                let w = [u * u * u, 3.0 * u * u * t, 3.0 * u * t * t, t * t * t];
                let dw = [
                    -3.0 * u * u,
                    3.0 * u * (u - 2.0 * t),
                    3.0 * t * (2.0 * u - t),
                    3.0 * t * t,
                ];
                let sum = |ws: [f64; 4], i: usize| -> f64 {
                    ws.iter()
                        .zip(p)
                        .map(|(w, q)| w * if i == 0 { q.0 } else { q.1 })
                        .sum()
                };
                ((sum(w, 0), sum(w, 1)), sum(dw, 0).hypot(sum(dw, 1)))
            }
        }
    }

    /// A wire of the bend's one edge, `mirror` as for [`Self::at`].
    fn wire(self, model: &mut Model, mirror: f64) -> Shape {
        let (curve, range): (ogeom_geom::Curve, (f64, f64)) = match self {
            Bend::Arc((x, y), r, (a, b)) => {
                let along = Direction::new(Vector::new(mirror, 0.0, 0.0), T).unwrap();
                let up = Direction::new(Vector::new(0.0, 0.0, mirror), T).unwrap();
                let frame = ogeom_math::Frame::new(Point::new(x, y, 0.0), up, along, T).unwrap();
                let circle = ogeom_math::Circle::new(frame, r, T).unwrap();
                (ogeom_geom::CircleCurve::new(circle).into(), (a, b))
            }
            Bend::Bezier(p) => (
                ogeom_geom::Curve::BSpline(
                    ogeom_geom::BSplineCurve::new(
                        ogeom_math::KnotVector::new(
                            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                            3,
                        )
                        .unwrap(),
                        p.iter().map(|&(x, y)| Point::new(x, y, 0.0)).collect(),
                        T,
                    )
                    .unwrap(),
                ),
                (0.0, 1.0),
            ),
        };
        let edge = ogeom_algo::make_edge(model, curve, range, T).unwrap().shape;
        make_wire(model, &[edge], T).unwrap().shape
    }
}

/// A curved edge of the profile inside a corner against an arc turning in
/// the arc's plane: its walls meet on a curve with no closed form, fitted
/// through each of its points' exact crossings, and both walls are cut on
/// that one edge. Each point's path is the closed form of the straight
/// profiles' test, integrated along the curve (split where it crosses the
/// corner's axis, where the path's length has a kink): arcs wholly inside
/// and crossing the axis twice, and Béziers wholly inside and across it,
/// mitred and round.
#[test]
fn a_curved_profile_edge_inside_a_corner_against_an_arc_sweeps_to_its_closed_form() {
    let bends = [
        Bend::Arc((0.5, 0.0), 1.0, (-FRAC_PI_2, FRAC_PI_2)),
        Bend::Arc((0.5, 0.0), 1.0, (FRAC_PI_2, 3.0 * FRAC_PI_2)),
        Bend::Bezier([(0.3, -1.0), (1.8, -0.6), (1.6, 0.7), (0.4, 1.0)]),
        Bend::Bezier([(-1.0, -1.0), (0.5, -1.2), (1.5, 0.5), (-0.5, 1.2)]),
    ];
    for (case, (name, pieces, flipped, run)) in inside_cases().iter().enumerate() {
        // A spline's walls are trimmed slowly: two spines take them.
        let bends = if case < 2 { &bends[..] } else { &bends[..2] };
        for &bend in bends {
            let (bend, mirror) = if *flipped {
                (bend.mirrored(), -1.0)
            } else {
                (bend, 1.0)
            };
            for corners in [PipeCorners::Mitre, PipeCorners::Round] {
                let round = corners == PipeCorners::Round;
                let mut model = Model::new();
                let wire = bend.wire(&mut model, mirror);
                let path = pieced(&mut model, pieces);
                let sheet = sweep(&mut model, &wire, &path, corners);
                let (area, boundary) = measured(&model, &sheet);
                // The kinks where the curve crosses `x = 0`, bisected.
                let x = |t: f64| bend.at(t, mirror).0.0;
                let mut breaks = vec![0.0];
                for k in 0..1024 {
                    let (mut a, mut b) = (f64::from(k) / 1024.0, f64::from(k + 1) / 1024.0);
                    if x(a) * x(b) < 0.0 {
                        for _ in 0..60 {
                            let m = f64::midpoint(a, b);
                            if x(a) * x(m) <= 0.0 {
                                b = m;
                            } else {
                                a = m;
                            }
                        }
                        breaks.push(f64::midpoint(a, b));
                    }
                }
                breaks.push(1.0);
                let along = |t: f64| {
                    let ((x, y), speed) = bend.at(t, mirror);
                    run(x, y, round) * speed
                };
                let speed = |t: f64| bend.at(t, mirror).1;
                let (mut want, mut length) = (0.0, 0.0);
                for w in breaks.windows(2) {
                    want += integral(along, w[0], w[1]);
                    length += integral(speed, w[0], w[1]);
                }
                let ends = [bend.at(0.0, mirror).0, bend.at(1.0, mirror).0];
                let around = 2.0 * length
                    + run(ends[0].0, ends[0].1, round)
                    + run(ends[1].0, ends[1].1, round);
                let case = format!("{name} {corners:?}");
                // Each wall is cut on a fitted curve held to the confusion
                // distance, which moves its area by at most that times the
                // cut's length: under the boundary's, for both walls.
                assert!(
                    (area - want).abs() <= 2.0 * T.confusion() * around,
                    "{case}: area {area} against {want}"
                );
                assert!(
                    close(boundary, around),
                    "{case}: boundary {boundary} against {around}"
                );
            }
        }
    }
}

/// A segment lying in the plane a corner against an arc turns in, reaching
/// inside it, sweeps a flat sheet cut where it crosses the corner's axis:
/// every face of it faces as the segment's prism down the first leg does,
/// and the other way for the segment reversed.
#[test]
fn an_open_sheet_inside_a_corner_against_an_arc_faces_as_its_prism() {
    let pieces = [
        Piece::Line(Point::new(0.0, 0.0, 10.0)),
        Piece::Arc {
            turned: Some(X),
            toward: Z,
            radius: 8.0,
            angle: FRAC_PI_2,
        },
    ];
    for reversed in [false, true] {
        for corners in [PipeCorners::Mitre, PipeCorners::Round] {
            let mut model = Model::new();
            let mut wire = polyline(&mut model, &[(-1.0, 0.0), (2.0, 0.0)]);
            if reversed {
                wire = wire.reversed();
            }
            let path = pieced(&mut model, &pieces);
            let sheet = sweep(&mut model, &wire, &path, corners);
            let prism = make_prism(&mut model, &wire, Vector::new(0.0, 0.0, 10.0), T)
                .unwrap()
                .shape;
            let face = explore(&model, &prism, Filter::OfType(ShapeType::Face)).unwrap();
            let (_, want) = ogeom_algo::face_normal(&model, &face[0], T).unwrap();
            let faces = explore(&model, &sheet, Filter::OfType(ShapeType::Face)).unwrap();
            assert!(faces.len() >= 4, "{reversed} {corners:?}");
            for face in &faces {
                let (at, n) = ogeom_algo::face_normal(&model, face, T).unwrap();
                assert!(
                    n.cross(want).magnitude() < 1e-9 && n.dot(want) > 0.0,
                    "{reversed} {corners:?}: {n:?} at {at:?} against {want:?}"
                );
            }
        }
    }
}

/// Inside a corner that turns an arc out of its plane a point's two paths
/// pass each other without meeting; a point `2` inside the turn onto an arc of radius 8 toward
/// `-z` meets it `atan(2 / √32) ≈ 0.34` rad round and `8 - √32 ≈ 2.34`
/// back from the corner, past an arc of `0.3` rad and a line of `2`: each
/// refused by name, mitred or round.
#[test]
fn an_open_profile_inside_a_corner_against_an_arc_it_cannot_join_is_refused_by_name() {
    let aside = vec![
        Piece::Line(Point::new(0.0, 0.0, 10.0)),
        Piece::Arc {
            turned: Some(X),
            toward: Y,
            radius: 8.0,
            angle: FRAC_PI_2,
        },
    ];
    for corners in [PipeCorners::Mitre, PipeCorners::Round] {
        let mut model = Model::new();
        let wire = polyline(&mut model, &L);
        let path = pieced(&mut model, &aside);
        let message = refusal(&mut model, &wire, &path, corners);
        assert!(
            message.contains(
                "an open wire reaches inside a corner that turns an arc out of its plane"
            ),
            "{corners:?}: {message}"
        );

        let deep = [(0.0, 0.0), (2.0, 0.0)];
        for (pieces, words) in [
            (
                vec![
                    Piece::Line(Point::new(0.0, 0.0, 10.0)),
                    Piece::Arc {
                        turned: Some(X),
                        toward: -Z,
                        radius: 8.0,
                        angle: 0.3,
                    },
                ],
                "meet past the arc's end",
            ),
            (
                vec![
                    Piece::Line(Point::new(0.0, 0.0, 2.0)),
                    Piece::Arc {
                        turned: Some(X),
                        toward: -Z,
                        radius: 8.0,
                        angle: FRAC_PI_2,
                    },
                ],
                "meet past the line leg's end",
            ),
        ] {
            let mut model = Model::new();
            let wire = polyline(&mut model, &deep);
            let path = pieced(&mut model, &pieces);
            let message = refusal(&mut model, &wire, &path, corners);
            assert!(message.contains(words), "{corners:?}: {message}");
        }
    }
}

/// A rectangle inside a corner against an arc turning in the arc's plane,
/// wholly or across the corner's axis, swept as a face, is bounded by its
/// two ends and the walls its two open halves sweep, within what the
/// face's skin is held to: the face ends each leg's walls inside the
/// corner where its points' paths cross, as the sheet does. A face wholly
/// inside has nothing to round, and is mitred only.
#[test]
fn two_open_halves_inside_a_corner_against_an_arc_bound_the_face_swept() {
    let pieces = [
        Piece::Line(Point::new(0.0, 0.0, 10.0)),
        Piece::Arc {
            turned: Some(X),
            toward: Z,
            radius: 8.0,
            angle: FRAC_PI_2,
        },
    ];
    for (rectangle, rounds) in [
        ([(0.5, -1.0), (2.0, -1.0), (2.0, 1.0), (0.5, 1.0)], false),
        ([(-1.0, -1.0), (2.0, -1.0), (2.0, 1.0), (-1.0, 1.0)], true),
    ] {
        let halves = [
            vec![rectangle[0], rectangle[1], rectangle[2]],
            vec![rectangle[2], rectangle[3], rectangle[0]],
        ];
        let (w, h) = (
            rectangle[1].0 - rectangle[0].0,
            rectangle[2].1 - rectangle[1].1,
        );
        let all = [PipeCorners::Mitre, PipeCorners::Round];
        for corners in &all[..if rounds { 2 } else { 1 }] {
            let mut model = Model::new();
            let path = pieced(&mut model, &pieces);
            let sheets: f64 = halves
                .iter()
                .map(|half| {
                    let wire = polyline(&mut model, half);
                    let sheet = sweep(&mut model, &wire, &path, *corners);
                    measured(&model, &sheet).0
                })
                .sum();
            let face = {
                let points: Vec<Point> = rectangle
                    .iter()
                    .map(|&(x, y)| Point::new(x, y, 0.0))
                    .collect();
                let wire = make_polygon(&mut model, &points, true, T).unwrap().shape;
                let plane = ogeom_math::Plane::new(
                    ogeom_math::Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap(),
                );
                ogeom_algo::make_face(
                    &mut model,
                    ogeom_geom::PlaneSurface::new(plane).into(),
                    &[wire],
                    T,
                )
                .unwrap()
                .shape
            };
            let solid = sweep(&mut model, &face, &path, *corners);
            let walls =
                surface_properties(&model, &solid, Deflection::with_chord(1e-4).unwrap(), T)
                    .unwrap()
                    .mass
                    - 2.0 * w * h;
            // The skin stands within its tolerance of the walls, whose
            // boundary is the rectangle's perimeter twice and the paths of
            // its four corners, each under 30 long.
            let bound = if *corners == PipeCorners::Mitre {
                1e-3 * (4.0 * (w + h) + 4.0 * 30.0)
            } else {
                walls * EXACT
            };
            assert!(
                (sheets - walls).abs() <= bound,
                "{rectangle:?} {corners:?}: halves {sheets} against the face's walls {walls}"
            );
        }
    }
}
