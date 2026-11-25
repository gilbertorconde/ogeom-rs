//! Largest-size readings confirmed against the build (ADR 0001, amendment
//! "largest-size readings confirmed by the build"): for each case the
//! reading's limit `L` builds valid at `L (1 - tol)` and is refused at
//! `L (1 + tol)`, `tol` the reading's own [`LIMIT_TOLERANCE`].
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{
    BlendLimit, BlendSizing, Chamfer, ChamferCorner, FilletCorner, LIMIT_TOLERANCE, RoundProfile,
    blend_limit,
};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A request and its options, owned.
enum Ask {
    Rounds(Vec<Shape>, FilletCorner),
    Distances(Shape, Shape, f64, f64),
    Chamfer(Shape, Chamfer, ChamferCorner),
}

impl Ask {
    fn sizing(&self) -> BlendSizing<'_> {
        match self {
            Self::Rounds(edges, corner) => BlendSizing::Rounds {
                edges,
                corner: *corner,
            },
            Self::Distances(edge, face, on_face, on_other) => BlendSizing::RoundDistances {
                edge,
                face,
                on_face: *on_face,
                on_other: *on_other,
                profile: RoundProfile::Elliptical,
            },
            Self::Chamfer(edge, spec, corner) => BlendSizing::Chamfer {
                edge,
                spec,
                corner: *corner,
            },
        }
    }

    /// The public call at `size`, on a copy of the model.
    fn build(&self, m: &Model, s: &Shape, size: f64) -> Result<(Model, Shape), String> {
        let mut m = m.clone();
        let built = match self {
            Self::Rounds(edges, corner) => {
                ogeom_fillet::fillet_edges_corner(&mut m, s, edges, size, *corner, T)
            }
            Self::Distances(edge, face, on_face, on_other) => ogeom_fillet::fillet_edge_distances(
                &mut m,
                s,
                edge,
                face,
                size * on_face,
                size * on_other,
                RoundProfile::Elliptical,
                T,
            ),
            Self::Chamfer(edge, spec, corner) => {
                let spec = match spec {
                    Chamfer::Symmetric(_) => Chamfer::Symmetric(size),
                    Chamfer::Distances {
                        face,
                        on_face,
                        on_other,
                    } => Chamfer::Distances {
                        face: face.clone(),
                        on_face: size * on_face,
                        on_other: size * on_other,
                    },
                    Chamfer::Angle { face, angle, .. } => Chamfer::Angle {
                        face: face.clone(),
                        distance: size,
                        angle: *angle,
                    },
                };
                ogeom_fillet::chamfer_edges_cornered(&mut m, s, &[(edge.clone(), spec)], *corner, T)
            }
        };
        built.map(|b| (m, b.shape)).map_err(|e| e.to_string())
    }
}

struct Case {
    name: &'static str,
    model: Model,
    solid: Shape,
    ask: Ask,
    up_to: f64,
}

/// The reading for `case`, held to its promise: valid at `L (1 - tol)`,
/// refused at `L (1 + tol)` (where the bound did not build); the refusal's
/// words returned.
fn confirmed(case: &Case) -> (BlendLimit, Option<String>) {
    let found = blend_limit(&case.model, &case.solid, &case.ask.sizing(), case.up_to, T)
        .unwrap_or_else(|e| panic!("{}: {e}", case.name));
    let l = found
        .size
        .unwrap_or_else(|| panic!("{}: no size builds: {found:?}", case.name));
    let (m, built) = case
        .ask
        .build(&case.model, &case.solid, l * (1.0 - LIMIT_TOLERANCE))
        .unwrap_or_else(|e| panic!("{}: refused below the limit {l}: {e}", case.name));
    let diagnosis = ogeom_algo::check(&m, &built, T).unwrap();
    assert!(
        diagnosis.is_valid(),
        "{}: invalid below the limit {l}: {:?}",
        case.name,
        diagnosis.problems
    );
    if found.refused.is_none() {
        assert_eq!(l, case.up_to, "{}", case.name);
        return (found, None);
    }
    let past = case
        .ask
        .build(&case.model, &case.solid, l * (1.0 + LIMIT_TOLERANCE));
    let refusal = past.err().unwrap_or_else(|| {
        panic!(
            "{}: builds a tolerance past the limit {l}: {found:?}",
            case.name
        )
    });
    (found, Some(refusal))
}

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// The edge of `s` passing nearest `p`, by 64 samples along each edge.
fn edge_near(m: &Model, s: &Shape, p: [f64; 3]) -> Shape {
    let q = Point::new(p[0], p[1], p[2]);
    let d = |e: &Shape| {
        let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
        (0..=64)
            .map(|i| {
                c.point_at(a + (b - a) * f64::from(i) / 64.0, T)
                    .unwrap()
                    .distance(q)
            })
            .fold(f64::INFINITY, f64::min)
    };
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| d(a).total_cmp(&d(b)))
        .unwrap()
}

/// The face of `s` whose bounds' centre is nearest `p`.
fn face_near(m: &Model, s: &Shape, p: [f64; 3]) -> Shape {
    let q = Point::new(p[0], p[1], p[2]);
    let d = |f: &Shape| {
        let b = ogeom_algo::shape_bounds(m, f, T).unwrap();
        let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
        Point::new(
            f64::midpoint(lo.x, hi.x),
            f64::midpoint(lo.y, hi.y),
            f64::midpoint(lo.z, hi.z),
        )
        .distance(q)
    };
    explore_unique(m, s, ShapeType::Face)
        .unwrap()
        .into_iter()
        .min_by(|a, b| d(a).total_cmp(&d(b)))
        .unwrap()
}

/// A prism of the (y, z) profile swept `length` along x.
fn prism(m: &mut Model, profile: &[(f64, f64)], length: f64) -> Shape {
    let points: Vec<Point> = profile
        .iter()
        .map(|&(y, z)| Point::new(0.0, y, z))
        .collect();
    let wire = ogeom_algo::make_polygon(m, &points, true, T).unwrap().shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::X),
        (-1000.0, 1000.0),
        (-1000.0, 1000.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let swept = ogeom_algo::make_prism(m, &face, Vector::new(length, 0.0, 0.0), T)
        .unwrap()
        .shape;
    ogeom_algo::baked_shape(m, &swept, T).unwrap().shape
}

fn boxed(dims: (f64, f64, f64), frame: Frame) -> (Model, Shape) {
    let mut m = Model::new();
    let s = ogeom_algo::make_box(&mut m, frame, dims, T).unwrap().shape;
    (m, s)
}

/// A frame leaning off every axis, as the stress suite places its parts.
fn oblique() -> Frame {
    let z = Direction::new(Vector::new(0.3, -0.5, 0.8), T).unwrap();
    let x = Direction::new(Vector::new(0.8, 0.6, 0.0).cross(z.vector()), T).unwrap();
    Frame::new(Point::new(1.3, -0.7, 2.1), z, x, T).unwrap()
}

fn all_edges(m: &Model, s: &Shape) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge).unwrap()
}

/// The edges along the sweep of a prism `length` long.
fn along(m: &Model, s: &Shape, length: f64) -> Vec<Shape> {
    all_edges(m, s)
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            let (p, q) = (c.point_at(a, T).unwrap(), c.point_at(b, T).unwrap());
            ((p.x - q.x).abs() - length).abs() < 1e-9
        })
        .collect()
}

fn plate() -> (Model, Shape) {
    boxed((20.0, 20.0, 2.0), at(0.0, 0.0, 0.0))
}

fn stepped() -> (Model, Shape) {
    let mut m = Model::new();
    let s = prism(
        &mut m,
        &[
            (0.0, 0.0),
            (30.0, 0.0),
            (30.0, 10.0),
            (15.0, 10.0),
            (15.0, 20.0),
            (0.0, 20.0),
        ],
        25.0,
    );
    (m, s)
}

fn cases() -> Vec<Case> {
    let mut cases = Vec::new();
    {
        let (m, s) = plate();
        let e = edge_near(&m, &s, [10.0, 0.0, 2.0]);
        cases.push(Case {
            name: "plate edge round",
            ask: Ask::Rounds(vec![e.clone()], FilletCorner::Setback),
            model: m.clone(),
            solid: s.clone(),
            up_to: 10.0,
        });
        cases.push(Case {
            name: "plate edge chamfer",
            ask: Ask::Chamfer(e.clone(), Chamfer::Symmetric(1.0), ChamferCorner::Facet),
            model: m.clone(),
            solid: s.clone(),
            up_to: 30.0,
        });
        let top = face_near(&m, &s, [10.0, 10.0, 2.0]);
        cases.push(Case {
            name: "plate edge chamfer distances",
            ask: Ask::Chamfer(
                e.clone(),
                Chamfer::Distances {
                    face: top.clone(),
                    on_face: 3.0,
                    on_other: 1.0,
                },
                ChamferCorner::Facet,
            ),
            model: m.clone(),
            solid: s.clone(),
            up_to: 30.0,
        });
        cases.push(Case {
            name: "plate edge distances round",
            ask: Ask::Distances(e.clone(), top, 3.0, 1.0),
            model: m.clone(),
            solid: s.clone(),
            up_to: 30.0,
        });
        let all = all_edges(&m, &s);
        cases.push(Case {
            name: "plate all round",
            ask: Ask::Rounds(all.clone(), FilletCorner::Setback),
            model: m.clone(),
            solid: s.clone(),
            up_to: 10.0,
        });
        cases.push(Case {
            name: "plate all round rolling",
            ask: Ask::Rounds(all, FilletCorner::Rolling),
            model: m.clone(),
            solid: s.clone(),
            up_to: 10.0,
        });
        let corner: Vec<Shape> = all_edges(&m, &s)
            .into_iter()
            .filter(|e| {
                let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
                [a, b].iter().any(|v| {
                    let p = m.node(v).unwrap().data().as_vertex().unwrap().point;
                    p.distance(Point::new(0.0, 0.0, 2.0)) < 1e-9
                })
            })
            .collect();
        cases.push(Case {
            name: "plate corner three round",
            ask: Ask::Rounds(corner, FilletCorner::Setback),
            model: m,
            solid: s,
            up_to: 10.0,
        });
    }
    {
        let (m, s) = boxed((3.1, 4.7, 2.3), oblique());
        let all = all_edges(&m, &s);
        cases.push(Case {
            name: "oblique box one edge round",
            ask: Ask::Rounds(vec![all[3].clone()], FilletCorner::Setback),
            model: m.clone(),
            solid: s.clone(),
            up_to: 10.0,
        });
        cases.push(Case {
            name: "oblique box all round",
            ask: Ask::Rounds(all.clone(), FilletCorner::Setback),
            model: m.clone(),
            solid: s.clone(),
            up_to: 10.0,
        });
        cases.push(Case {
            name: "oblique box one edge chamfer",
            ask: Ask::Chamfer(
                all[5].clone(),
                Chamfer::Symmetric(1.0),
                ChamferCorner::Facet,
            ),
            model: m.clone(),
            solid: s.clone(),
            up_to: 10.0,
        });
    }
    {
        let mut m = Model::new();
        let s = ogeom_algo::make_cylinder(&mut m, at(0.0, 0.0, 0.0), 2.5, 3.0, T)
            .unwrap()
            .shape;
        let rim = edge_near(&m, &s, [0.0, 2.5, 3.0]);
        let low = edge_near(&m, &s, [0.0, 2.5, 0.0]);
        cases.push(Case {
            name: "drum rim round",
            ask: Ask::Rounds(vec![rim.clone()], FilletCorner::Setback),
            model: m.clone(),
            solid: s.clone(),
            up_to: 10.0,
        });
        cases.push(Case {
            name: "drum both rims round",
            ask: Ask::Rounds(vec![rim.clone(), low], FilletCorner::Setback),
            model: m.clone(),
            solid: s.clone(),
            up_to: 10.0,
        });
        cases.push(Case {
            name: "drum rim chamfer",
            ask: Ask::Chamfer(rim, Chamfer::Symmetric(1.0), ChamferCorner::Facet),
            model: m,
            solid: s,
            up_to: 10.0,
        });
    }
    {
        let mut m = Model::new();
        let block = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (6.0, 6.0, 3.0), T)
            .unwrap()
            .shape;
        let drill = ogeom_algo::make_cylinder(&mut m, at(3.0, 3.0, -1.0), 1.2, 5.0, T)
            .unwrap()
            .shape;
        let s = ogeom_bool::cut(&mut m, &block, &drill, T).unwrap().shape;
        let rim = edge_near(&m, &s, [3.0, 4.2, 3.0]);
        cases.push(Case {
            name: "bore rim round",
            ask: Ask::Rounds(vec![rim.clone()], FilletCorner::Setback),
            model: m.clone(),
            solid: s.clone(),
            up_to: 10.0,
        });
        cases.push(Case {
            name: "bore rim chamfer",
            ask: Ask::Chamfer(rim, Chamfer::Symmetric(1.0), ChamferCorner::Facet),
            model: m,
            solid: s,
            up_to: 10.0,
        });
    }
    {
        let mut m = Model::new();
        let s = ogeom_algo::make_cone(&mut m, Frame::WORLD, 6.0, 5.0, 2.0, T)
            .unwrap()
            .shape;
        let rim = edge_near(&m, &s, [0.0, 5.0, 2.0]);
        cases.push(Case {
            name: "cone rim round",
            ask: Ask::Rounds(vec![rim], FilletCorner::Setback),
            model: m,
            solid: s,
            up_to: 10.0,
        });
    }
    {
        let (m, s) = stepped();
        let riser_top = edge_near(&m, &s, [12.0, 15.0, 20.0]);
        cases.push(Case {
            name: "step riser top chamfer",
            ask: Ask::Chamfer(
                riser_top.clone(),
                Chamfer::Symmetric(1.0),
                ChamferCorner::Facet,
            ),
            model: m.clone(),
            solid: s.clone(),
            up_to: 30.0,
        });
        cases.push(Case {
            name: "step riser top round",
            ask: Ask::Rounds(vec![riser_top], FilletCorner::Setback),
            model: m.clone(),
            solid: s.clone(),
            up_to: 30.0,
        });
        // The riser's foot (concave) with the two convex side edges at its
        // end: a mixed corner.
        let foot = edge_near(&m, &s, [12.0, 15.0, 10.0]);
        let upright = edge_near(&m, &s, [0.0, 15.0, 15.0]);
        let floor_side = edge_near(&m, &s, [0.0, 22.0, 10.0]);
        cases.push(Case {
            name: "step mixed corner round",
            ask: Ask::Rounds(
                vec![foot.clone(), upright.clone(), floor_side.clone()],
                FilletCorner::Setback,
            ),
            model: m.clone(),
            solid: s.clone(),
            up_to: 30.0,
        });
        cases.push(Case {
            name: "step foot and upright round",
            ask: Ask::Rounds(vec![foot.clone(), upright], FilletCorner::Setback),
            model: m.clone(),
            solid: s.clone(),
            up_to: 30.0,
        });
        cases.push(Case {
            name: "step foot and floor side round",
            ask: Ask::Rounds(vec![foot, floor_side], FilletCorner::Setback),
            model: m.clone(),
            solid: s.clone(),
            up_to: 30.0,
        });
        cases.push(Case {
            name: "stepped block all round",
            ask: Ask::Rounds(all_edges(&m, &s), FilletCorner::Setback),
            model: m,
            solid: s,
            up_to: 30.0,
        });
    }
    {
        // Two top edges at a plate's corner, the upright between them
        // left sharp.
        let (m, s) = plate();
        let a = edge_near(&m, &s, [10.0, 0.0, 2.0]);
        let b = edge_near(&m, &s, [0.0, 10.0, 2.0]);
        cases.push(Case {
            name: "plate two of three round",
            ask: Ask::Rounds(vec![a.clone(), b.clone()], FilletCorner::Setback),
            model: m.clone(),
            solid: s.clone(),
            up_to: 10.0,
        });
        let c = edge_near(&m, &s, [0.0, 0.0, 1.0]);
        cases.push(Case {
            name: "plate corner rolling",
            ask: Ask::Rounds(vec![a, b, c], FilletCorner::Rolling),
            model: m,
            solid: s,
            up_to: 10.0,
        });
    }
    for (name, every, profile) in [
        (
            "L channel along round",
            "L channel every edge round",
            vec![
                (0.0, 0.0),
                (20.0, 0.0),
                (20.0, 5.0),
                (5.0, 5.0),
                (5.0, 20.0),
                (0.0, 20.0),
            ],
        ),
        (
            "U channel along round",
            "U channel every edge round",
            vec![
                (0.0, 0.0),
                (30.0, 0.0),
                (30.0, 20.0),
                (25.0, 20.0),
                (25.0, 5.0),
                (5.0, 5.0),
                (5.0, 20.0),
                (0.0, 20.0),
            ],
        ),
    ] {
        let mut m = Model::new();
        let s = prism(&mut m, &profile, 30.0);
        let edges = along(&m, &s, 30.0);
        cases.push(Case {
            name,
            ask: Ask::Rounds(edges, FilletCorner::Setback),
            model: m.clone(),
            solid: s.clone(),
            up_to: 30.0,
        });
        cases.push(Case {
            name: every,
            ask: Ask::Rounds(all_edges(&m, &s), FilletCorner::Setback),
            model: m,
            solid: s,
            up_to: 30.0,
        });
    }
    {
        // A boss standing on a block near its edge: the boss's foot round
        // reaching the block's top edge.
        let mut m = Model::new();
        let block = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (20.0, 20.0, 5.0), T)
            .unwrap()
            .shape;
        let boss = ogeom_algo::make_cylinder(&mut m, at(4.0, 10.0, 5.0), 2.0, 4.0, T)
            .unwrap()
            .shape;
        let s = ogeom_bool::fuse(&mut m, &block, &boss, T).unwrap().shape;
        let foot = edge_near(&m, &s, [4.0, 12.0, 5.0]);
        cases.push(Case {
            name: "boss foot round",
            ask: Ask::Rounds(vec![foot], FilletCorner::Setback),
            model: m,
            solid: s,
            up_to: 10.0,
        });
    }
    {
        // A cube with one corner chopped off square to its diagonal: three
        // oblique trihedral corners round the triangle.
        let mut m = Model::new();
        let cube = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let n = Direction::new(Vector::new(1.0, 1.0, 1.0), T).unwrap();
        let fx = Direction::new(Vector::new(1.0, -1.0, 0.0), T).unwrap();
        let fy = Direction::new(n.vector().cross(fx.vector()), T).unwrap();
        let p = Point::new(7.0, 10.0, 10.0);
        let origin = p - fx.vector() * 30.0 - fy.vector() * 30.0;
        let knife = ogeom_algo::make_box(
            &mut m,
            Frame::new(origin, n, fx, T).unwrap(),
            (60.0, 60.0, 30.0),
            T,
        )
        .unwrap()
        .shape;
        let s = ogeom_bool::cut(&mut m, &cube, &knife, T).unwrap().shape;
        let near: Vec<Shape> = all_edges(&m, &s)
            .into_iter()
            .filter(|e| {
                let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
                [a, b].iter().any(|v| {
                    let q = m.node(v).unwrap().data().as_vertex().unwrap().point;
                    q.x + q.y + q.z > 26.9
                })
            })
            .collect();
        cases.push(Case {
            name: "chopped cube corner round",
            ask: Ask::Rounds(near, FilletCorner::Setback),
            model: m.clone(),
            solid: s.clone(),
            up_to: 10.0,
        });
        cases.push(Case {
            name: "chopped cube all round",
            ask: Ask::Rounds(all_edges(&m, &s), FilletCorner::Setback),
            model: m,
            solid: s,
            up_to: 10.0,
        });
    }
    cases
}

/// The cases on curved faces, and the corpus's prismatic parts, each
/// taking a second or more.
const HEAVY: [&str; 15] = [
    "drum rim round",
    "drum both rims round",
    "drum rim chamfer",
    "bore rim round",
    "bore rim chamfer",
    "cone rim round",
    "step mixed corner round",
    "step foot and upright round",
    "stepped block all round",
    "L channel along round",
    "L channel every edge round",
    "U channel along round",
    "U channel every edge round",
    "boss foot round",
    "chopped cube all round",
];

/// Every case of one weight confirmed, by name.
fn swept(heavy: bool) -> Vec<(&'static str, BlendLimit, Option<String>)> {
    cases()
        .iter()
        .filter(|c| HEAVY.contains(&c.name) == heavy)
        .map(|c| {
            let (found, refusal) = confirmed(c);
            (c.name, found, refusal)
        })
        .collect()
}

fn find<'a>(
    found: &'a [(&'static str, BlendLimit, Option<String>)],
    name: &str,
) -> (f64, &'a BlendLimit, &'a str) {
    let (_, limit, refusal) = found.iter().find(|(n, ..)| *n == name).unwrap();
    (limit.size.unwrap(), limit, refusal.as_deref().unwrap_or(""))
}

/// `size` stands within the reading's tolerance below `exact`, or a
/// confusion distance under it where the seat reading stands short.
fn near(name: &str, size: f64, exact: f64) {
    assert!(
        size <= exact * (1.0 + LIMIT_TOLERANCE) && size >= exact * (1.0 - 2.0 * LIMIT_TOLERANCE),
        "{name}: limit {size} against {exact}"
    );
}

/// The planar cases. Refused for each reason the planar engine has:
/// - setback: a 2 mm plate's edge rounds to its thickness;
/// - a corner that does not close: a plate's corner of three rounds
///   rolling, refused as soon as a round reaches the floor beyond its
///   ends, short of the setback; and a cube's corner chopped off square to
///   its diagonal, whose three oblique corners stop the rounds at 2.366,
///   short of what one of their edges rounds to alone;
/// - a concave floor: a step's riser top bevels at most the riser's 10 mm
///   height, past which the bevel would cut below the floor;
/// - self-intersection: a plate's edge bevels, rolling over the side's
///   foot as a notch, until the bevel's setbacks cross on the end faces,
///   where its plane passes the plate's far bottom edge: 20 + 2 for the
///   symmetric bevel, 20/3 + 2 for the scale of 3 along the top and 1 down
///   the side.
///
/// Three rounds at the plate's corner build up to the plate's thickness,
/// through the one-boolean fallback, whose polytope's rounds reach the
/// floor beyond their ends as the plate's do.
#[test]
fn planar_limits_build_below_and_are_refused_above() {
    let found = swept(false);
    let (size, limit, why) = find(&found, "plate edge round");
    near("plate edge round", size, 2.0);
    assert_eq!(
        limit.probes,
        2 + 2 * ogeom_fillet::LIMIT_SAMPLES,
        "{limit:?}"
    );
    assert!(why.contains("does not reach so far"), "{why}");

    let (size, _, why) = find(&found, "plate corner rolling");
    let (setback, ..) = find(&found, "plate corner three round");
    near("plate corner three round", setback, 2.0);
    assert!(size < setback, "{size} against {setback}");
    assert!(why.contains("rolling corner"), "{why}");
    let (size, limit, why) = find(&found, "chopped cube corner round");
    assert!(
        limit.seat.unwrap() > 1.4 * size,
        "the corner stops short of the seat: {limit:?}"
    );
    assert!(why.contains("builds no round of this kind"), "{why}");

    let (size, _, why) = find(&found, "step riser top chamfer");
    near("step riser top chamfer", size, 10.0);
    assert!(why.contains("concave edge"), "{why}");

    let (size, limit, why) = find(&found, "plate edge chamfer");
    near("plate edge chamfer", size, 22.0);
    assert_eq!(limit.seat, Some(2.0), "{limit:?}");
    assert!(why.contains("setbacks on it cross"), "{why}");
    let (size, _, why) = find(&found, "plate edge chamfer distances");
    near("plate edge chamfer distances", size, 20.0 / 3.0 + 2.0);
    assert!(why.contains("setbacks on it cross"), "{why}");

    for (name, limit, _) in &found {
        assert_eq!(limit.builds_again, None, "{name}: {limit:?}");
    }
}

/// The curved cases and the corpus's parts. Refused for the reasons the
/// planar cases do not reach:
/// - validity: both rims of a drum 3 high round until their bands meet
///   halfway up the wall, where the result fails its check;
/// - overflow: a boss of radius 2 standing 2 in from a block's edge has
///   its foot rounded until the round's contact runs over that edge, the
///   overflow the engine does not build.
#[test]
#[ignore = "heavy"]
fn curved_and_corpus_limits_build_below_and_are_refused_above() {
    let found = swept(true);
    let (size, _, why) = find(&found, "drum both rims round");
    near("drum both rims round", size, 1.5);
    assert!(why.contains("fails its check"), "{why}");

    let (size, _, why) = find(&found, "boss foot round");
    assert!((2.0..2.01).contains(&size), "boss foot limit {size}");
    assert!(why.contains("builds no round of this kind"), "{why}");

    // The sweep edges of the L and U channels, 5 thick: the bands either
    // side of a 5 mm wall meet at 2.5.
    for name in ["L channel along round", "U channel along round"] {
        near(name, find(&found, name).0, 2.5);
    }
    for (name, limit, _) in &found {
        assert_eq!(limit.builds_again, None, "{name}: {limit:?}");
    }
}

/// One edge of the chopped cube's corner triangle rounds alone past where
/// the corner stops the three: the refusal at the corner is the corner's.
#[test]
fn the_chopped_corner_s_edge_alone_rounds_past_the_corner() {
    let all = cases();
    let corner = all
        .iter()
        .find(|c| c.name == "chopped cube corner round")
        .unwrap();
    let (together, _) = confirmed(corner);
    let Ask::Rounds(edges, _) = &corner.ask else {
        unreachable!()
    };
    let alone = edges
        .iter()
        .map(|e| {
            blend_limit(
                &corner.model,
                &corner.solid,
                &BlendSizing::Rounds {
                    edges: std::slice::from_ref(e),
                    corner: FilletCorner::Setback,
                },
                10.0,
                T,
            )
            .unwrap()
            .size
            .unwrap()
        })
        .fold(f64::INFINITY, f64::min);
    assert!(
        alone > together.size.unwrap() * 1.05,
        "alone {alone}, together {together:?}"
    );
}
