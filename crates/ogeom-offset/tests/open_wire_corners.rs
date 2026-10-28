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

/// A point inside a corner against an arc crosses the plane halving it on
/// the legs themselves, where the line's path and the arc's do not meet
/// it at one point: refused by name, mitred or round, before or after the
/// arc.
#[test]
fn an_open_profile_inside_a_corner_against_an_arc_is_refused_by_name() {
    let pieces = [
        vec![
            Piece::Line(Point::new(0.0, 0.0, 10.0)),
            Piece::Arc {
                turned: Some(X),
                toward: Z,
                radius: 8.0,
                angle: FRAC_PI_2,
            },
        ],
        vec![
            Piece::Arc {
                turned: None,
                toward: X,
                radius: 8.0,
                angle: FRAC_PI_2,
            },
            Piece::Line(Point::new(8.0, 0.0, 18.0)),
        ],
    ];
    for pieces in &pieces {
        for corners in [PipeCorners::Mitre, PipeCorners::Round] {
            let mut model = Model::new();
            let wire = polyline(&mut model, &L);
            let path = pieced(&mut model, pieces);
            let message = refusal(&mut model, &wire, &path, corners);
            assert!(
                message.contains("an open wire reaches inside a corner against an arc"),
                "{corners:?}: {message}"
            );
        }
    }
}

/// The binormal square to a plane spine, and the Frenet frame of two arcs
/// whose normals the corner between them turns onto each other (a quarter
/// arc toward `+x`, then sharply up `+z` into one toward `+x` again), are
/// the rotation-minimizing frame turned by a fixed angle: the open sheet
/// is the one that frame sweeps. A binormal leaning off the spine's plane
/// sets the next leg's section other than the last turned about the
/// corner, as does the Frenet frame of a second arc bending toward `-x`,
/// whose normal the corner turns the other way, and a straight run between
/// corners has no Frenet frame: each refused by name.
#[test]
fn a_law_carrying_the_section_as_the_rotation_minimizing_frame_sweeps_its_sheet() {
    let l_spine = [
        Point::ORIGIN,
        Point::new(0.0, 0.0, 10.0),
        Point::new(10.0, 0.0, 10.0),
    ];
    let arcs = |toward: Vector| {
        [
            Piece::Arc {
                turned: None,
                toward: X,
                radius: 8.0,
                angle: FRAC_PI_2,
            },
            Piece::Arc {
                turned: Some(Z),
                toward,
                radius: 5.0,
                angle: FRAC_PI_2,
            },
        ]
    };
    let leaning = PipeLaw::Binormal(Direction::new(Vector::new(1.0, 1.0, 0.0), T).unwrap());
    let level = PipeLaw::Binormal(Direction::Y);
    for corners in [PipeCorners::Mitre, PipeCorners::Round] {
        let mut model = Model::new();
        let wire = polyline(&mut model, &L);
        let path = spine(&mut model, &l_spine);
        let carried = sweep(&mut model, &wire, &path, corners);
        let want = measured(&model, &carried);
        let by_law = make_pipe_shell_with(&mut model, &wire, &path, &level, corners, 1e-3, T)
            .unwrap_or_else(|e| panic!("binormal {corners:?}: {e}"))
            .shape;
        let got = measured(&model, &by_law);
        assert!(
            close(got.0, want.0) && close(got.1, want.1),
            "binormal {corners:?}: {got:?} against {want:?}"
        );
        let message =
            match make_pipe_shell_with(&mut model, &wire, &path, &leaning, corners, 1e-3, T) {
                Ok(_) => panic!("leaning binormal {corners:?}: built"),
                Err(e) => e.to_string(),
            };
        assert!(
            message.contains("the law turns an open wire's section other than"),
            "leaning binormal {corners:?}: {message}"
        );
    }
    let mut model = Model::new();
    let wire = polyline(&mut model, &MIRRORED);
    let path = pieced(&mut model, &arcs(X));
    let carried = sweep(&mut model, &wire, &path, PipeCorners::Round);
    let want = measured(&model, &carried);
    let by_law = make_pipe_shell_with(
        &mut model,
        &wire,
        &path,
        &PipeLaw::Frenet,
        PipeCorners::Round,
        1e-3,
        T,
    )
    .unwrap()
    .shape;
    let got = measured(&model, &by_law);
    assert!(
        close(got.0, want.0) && close(got.1, want.1),
        "Frenet: {got:?} against {want:?}"
    );
    for (path, says) in [
        (
            pieced(&mut model, &arcs(-X)),
            "the law turns an open wire's section other than",
        ),
        (
            spine(&mut model, &l_spine),
            "the Frenet frame is not defined along a straight run",
        ),
    ] {
        let message = match make_pipe_shell_with(
            &mut model,
            &wire,
            &path,
            &PipeLaw::Frenet,
            PipeCorners::Round,
            1e-3,
            T,
        ) {
            Ok(_) => panic!("{says}: built"),
            Err(e) => e.to_string(),
        };
        assert!(message.contains(says), "{message}");
    }
}

/// A spine of a line and an arc placed by a location, on the wire or on
/// each edge under it, with the profile moved alike, sweeps the closed
/// form of the same shapes unmoved: up `+z` 10, then sharply `+x` into a
/// quarter arc of radius 8 toward `+y`, the corner `α = -x` and the arc's
/// axis `8 - y` off.
#[test]
fn an_open_sweep_against_an_arc_down_a_placed_spine_keeps_its_closed_form() {
    let motion = Transform::translation(Vector::new(3.0, -7.0, 2.5))
        * Transform::rotation(
            Axis {
                location: Point::new(1.0, 2.0, 3.0),
                direction: Direction::new(Vector::new(1.0, 2.0, -0.5), T).unwrap(),
            },
            1.1,
        );
    let pieces = [
        Piece::Line(Point::new(0.0, 0.0, 10.0)),
        Piece::Arc {
            turned: Some(X),
            toward: Y,
            radius: 8.0,
            angle: FRAC_PI_2,
        },
    ];
    for corners in [PipeCorners::Mitre, PipeCorners::Round] {
        let gain = |alpha: f64| {
            if corners == PipeCorners::Round {
                FRAC_PI_2 * alpha
            } else {
                2.0 * alpha
            }
        };
        let (want, around) =
            linear_closed_form(&OUTSIDE, |x, y| 10.0 + gain(-x) + FRAC_PI_2 * (8.0 - y));
        for on_edges in [false, true] {
            let mut model = Model::new();
            let profile: Vec<Point> = OUTSIDE
                .iter()
                .map(|&(x, y)| motion.apply(Point::new(x, y, 0.0)))
                .collect();
            let wire = make_polygon(&mut model, &profile, false, T).unwrap().shape;
            let placed = Location::of(model.add_datum(motion));
            let path = pieced(&mut model, &pieces);
            let path = if on_edges {
                let edges: Vec<Shape> = model
                    .ordered_children_of(&path)
                    .unwrap()
                    .iter()
                    .map(|e| e.moved(&placed))
                    .collect();
                make_wire(&mut model, &edges, T).unwrap().shape
            } else {
                path.moved(&placed)
            };
            let sheet = sweep(&mut model, &wire, &path, corners);
            let (area, boundary) = measured(&model, &sheet);
            let case = format!("{corners:?} on edges {on_edges}");
            assert!(close(area, want), "{case}: {area} against {want}");
            assert!(
                close(boundary, around),
                "{case}: {boundary} against {around}"
            );
        }
    }
}
