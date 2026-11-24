//! The largest measure several chamfers take at once: every spelling,
//! against closed forms where the refusal just past them is the overlap or
//! the setback, builds at the answer under both corners and is refused just
//! past it. Where the request builds past its seat, the seat is still the
//! closed form and the answer is past it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use std::f64::consts::{FRAC_PI_3, FRAC_PI_4, FRAC_PI_6};

use ogeom_core::Tolerances;
use ogeom_fillet::{
    BlendSizing, Chamfer, ChamferCorner, LIMIT_TOLERANCE, Overflow, blend_limit,
    chamfer_edges_cornered, largest_chamfer_distance_all,
};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const CORNERS: [ChamferCorner; 2] = [ChamferCorner::Facet, ChamferCorner::Point];

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

/// The face of `s` whose bounds' centre lies nearest `p`.
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

/// `specs` with every measure set to `measure`, as the reading reads it.
fn at_measure(specs: &[(Shape, Chamfer)], measure: f64) -> Vec<(Shape, Chamfer)> {
    specs
        .iter()
        .map(|(edge, spec)| {
            let scaled = match spec {
                Chamfer::Symmetric(_) => Chamfer::Symmetric(measure),
                Chamfer::Distances {
                    face,
                    on_face,
                    on_other,
                } => Chamfer::Distances {
                    face: face.clone(),
                    on_face: on_face * measure,
                    on_other: on_other * measure,
                },
                Chamfer::Angle { face, angle, .. } => Chamfer::Angle {
                    face: face.clone(),
                    distance: measure,
                    angle: *angle,
                },
            };
            (edge.clone(), scaled)
        })
        .collect()
}

/// The reading of `specs` lies within `1e-6 · scale` under `oracle`, and
/// the request builds a valid solid at the reading and at
/// `oracle - 1e-6 · scale` under both corners, with `faces` faces where
/// given.
fn reads_and_builds(
    m: &Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    oracle: f64,
    scale: f64,
    faces: Option<usize>,
) -> f64 {
    let got = largest_chamfer_distance_all(m, solid, specs, 100.0, T)
        .unwrap()
        .unwrap();
    assert!(
        (oracle - 1e-6 * scale..=oracle).contains(&got),
        "largest measure {got} against {oracle}"
    );
    for measure in [got, oracle - 1e-6 * scale] {
        for corner in CORNERS {
            let mut built_in = m.clone();
            let built = chamfer_edges_cornered(
                &mut built_in,
                solid,
                &at_measure(specs, measure),
                corner,
                T,
            )
            .unwrap_or_else(|e| panic!("{corner:?} at {measure}, read {got}: {e}"));
            let diagnosis = ogeom_algo::check(&built_in, &built.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{corner:?}: {:?}", diagnosis.problems);
            if let Some(faces) = faces {
                let n = explore_unique(&built_in, &built.shape, ShapeType::Face)
                    .unwrap()
                    .len();
                assert_eq!(n, faces, "{corner:?} at {measure}");
            }
        }
    }
    got
}

/// The request at `measure` is refused under both corners.
fn refused_at(m: &Model, solid: &Shape, specs: &[(Shape, Chamfer)], measure: f64) {
    for corner in CORNERS {
        let mut past = m.clone();
        assert!(
            chamfer_edges_cornered(&mut past, solid, &at_measure(specs, measure), corner, T)
                .is_err(),
            "{corner:?} at {measure}, past the largest measure, builds"
        );
    }
}

/// A bar 20 long, 2 wide and 10 high, and its two long top edges.
fn bar(m: &mut Model) -> (Shape, [Shape; 2], Shape) {
    let bar = ogeom_algo::make_box(m, at(0.0, 0.0, 0.0), (20.0, 2.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = [
        edge_near(m, &bar, [10.0, 0.0, 10.0]),
        edge_near(m, &bar, [10.0, 2.0, 10.0]),
    ];
    let top = face_near(m, &bar, [10.0, 1.0, 10.0]);
    (bar, edges, top)
}

/// The bevels of a narrow face's two edges meet across it at half its
/// width over their legs on it, every spelling: 1 symmetric, a scale of
/// 1/2 for 2 along the top and 1 down the side, and 1 at an angle from
/// the top. That is the seat: there the top survives as a sliver (8
/// faces), and just past it the two convex bevels meet along their ridge
/// and the top is gone (7 faces), which is built, not refused. So the
/// answer stands past the seat, where the build itself is refused.
#[test]
fn a_narrow_face_s_two_edges_chamfer_past_half_its_width() {
    let mut m = Model::new();
    let (bar, edges, top) = bar(&mut m);
    let spellings: [(Chamfer, f64); 3] = [
        (Chamfer::Symmetric(0.3), 1.0),
        (
            Chamfer::Distances {
                face: top.clone(),
                on_face: 2.0,
                on_other: 1.0,
            },
            0.5,
        ),
        (
            Chamfer::Angle {
                face: top.clone(),
                distance: 0.3,
                angle: FRAC_PI_3,
            },
            1.0,
        ),
    ];
    for (spec, oracle) in spellings {
        let specs: Vec<_> = edges.iter().map(|e| (e.clone(), spec.clone())).collect();
        let sizing = BlendSizing::Chamfers {
            specs: &specs,
            corner: ChamferCorner::Facet,
            overflow: Overflow::default(),
        };
        let found = blend_limit(&m, &bar, &sizing, 100.0, T).unwrap();
        let seat = found.seat.unwrap();
        assert!(
            (oracle - 1e-6 * 20.0..=oracle).contains(&seat),
            "{spec:?}: seat {found:?} against {oracle}"
        );
        for corner in CORNERS {
            for (measure, faces) in [(seat, 8), (seat * 1.01, 7)] {
                let mut built_in = m.clone();
                let built = chamfer_edges_cornered(
                    &mut built_in,
                    &bar,
                    &at_measure(&specs, measure),
                    corner,
                    T,
                )
                .unwrap();
                let n = explore_unique(&built_in, &built.shape, ShapeType::Face)
                    .unwrap()
                    .len();
                assert_eq!(n, faces, "{spec:?} {corner:?} at {measure}");
            }
        }
        let got = largest_chamfer_distance_all(&m, &bar, &specs, 100.0, T)
            .unwrap()
            .unwrap();
        assert_eq!(Some(got), found.size);
        assert!(got > seat * 1.01, "{spec:?}: {found:?}");
        let mut at = m.clone();
        let built = chamfer_edges_cornered(
            &mut at,
            &bar,
            &at_measure(&specs, got),
            ChamferCorner::Facet,
            T,
        )
        .unwrap();
        assert!(ogeom_algo::check(&at, &built.shape, T).unwrap().is_valid());
        let mut past = m.clone();
        assert!(
            chamfer_edges_cornered(
                &mut past,
                &bar,
                &at_measure(&specs, got * (1.0 + LIMIT_TOLERANCE)),
                ChamferCorner::Facet,
                T,
            )
            .is_err(),
            "{spec:?}: builds past {found:?}"
        );
    }
}

/// A block 20 x 10 x 10 with a slot 2 wide and 5 deep along it, and the
/// slot's two floor edges.
fn slot(m: &mut Model) -> (Shape, [Shape; 2], Shape) {
    let block = ogeom_algo::make_box(m, at(0.0, 0.0, 0.0), (20.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let tool = ogeom_algo::make_box(m, at(-1.0, 4.0, 5.0), (22.0, 2.0, 6.0), T)
        .unwrap()
        .shape;
    let slot = ogeom_bool::cut(m, &block, &tool, T).unwrap().shape;
    let edges = [
        edge_near(m, &slot, [10.0, 4.0, 5.0]),
        edge_near(m, &slot, [10.0, 6.0, 5.0]),
    ];
    let floor = face_near(m, &slot, [10.0, 5.0, 5.0]);
    (slot, edges, floor)
}

/// A slot of width 2: the concave bevels of its floor edges meet across
/// the floor at w/2 = 1 symmetric, a scale of 1/2 for 2 along the floor
/// and 1 up the wall, and 1 at an angle from the floor; just past it the
/// floor's setbacks cross and the request is refused.
#[test]
fn a_slot_s_floor_edges_chamfer_at_most_half_its_width() {
    let mut m = Model::new();
    let (slot, edges, floor) = slot(&mut m);
    let spellings: [(Chamfer, f64); 3] = [
        (Chamfer::Symmetric(0.3), 1.0),
        (
            Chamfer::Distances {
                face: floor.clone(),
                on_face: 2.0,
                on_other: 1.0,
            },
            0.5,
        ),
        (
            Chamfer::Angle {
                face: floor.clone(),
                distance: 0.3,
                angle: FRAC_PI_6,
            },
            1.0,
        ),
    ];
    for (spec, oracle) in spellings {
        let specs: Vec<_> = edges.iter().map(|e| (e.clone(), spec.clone())).collect();
        let got = reads_and_builds(&m, &slot, &specs, oracle, 20.0, None);
        refused_at(&m, &slot, &specs, got * 1.01);
    }
}

/// A post 10 x 10, 6 high, standing on a base 30 x 30 x 5, and the three
/// edges at one top corner of the post: along x, along y and down.
fn post(m: &mut Model) -> (Shape, [Shape; 3]) {
    let base = ogeom_algo::make_box(m, at(0.0, 0.0, 0.0), (30.0, 30.0, 5.0), T)
        .unwrap()
        .shape;
    let post = ogeom_algo::make_box(m, at(10.0, 10.0, 5.0), (10.0, 10.0, 6.0), T)
        .unwrap()
        .shape;
    let part = ogeom_bool::fuse(m, &base, &post, T).unwrap().shape;
    let edges = [
        edge_near(m, &part, [15.0, 10.0, 11.0]),
        edge_near(m, &part, [10.0, 15.0, 11.0]),
        edge_near(m, &part, [10.0, 10.0, 8.0]),
    ];
    (part, edges)
}

/// Three edges meeting at a corner are corner mates, not held apart, so
/// the limit is the single edges': where a top edge's leg down the post
/// reaches the floor 6 below. Symmetric, that is 6; with 1 along the top
/// (or the post's x = 10 side, for the upright edge) and 2 along the
/// other, a scale of 3; at 60 degrees from the top (45 from the side for
/// the upright edge), 6 / tan 60 = 2√3. Just past it the leg runs below
/// the floor and is refused, under both corners.
#[test]
fn three_edges_at_a_corner_chamfer_as_far_as_each_alone() {
    let mut m = Model::new();
    let (part, edges) = post(&mut m);
    let top = face_near(&m, &part, [15.0, 15.0, 11.0]);
    let side = face_near(&m, &part, [10.0, 15.0, 8.0]);
    let named = [top.clone(), top, side];
    let symmetric: Vec<_> = edges
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(0.3)))
        .collect();
    let distances: Vec<_> = edges
        .iter()
        .zip(&named)
        .map(|(e, f)| {
            (
                e.clone(),
                Chamfer::Distances {
                    face: f.clone(),
                    on_face: 1.0,
                    on_other: 2.0,
                },
            )
        })
        .collect();
    let angles: Vec<_> = edges
        .iter()
        .zip(&named)
        .zip([FRAC_PI_3, FRAC_PI_3, FRAC_PI_4])
        .map(|((e, f), angle)| {
            (
                e.clone(),
                Chamfer::Angle {
                    face: f.clone(),
                    distance: 0.3,
                    angle,
                },
            )
        })
        .collect();
    for (specs, oracle) in [
        (symmetric, 6.0),
        (distances, 3.0),
        (angles, 2.0 * 3.0_f64.sqrt()),
    ] {
        let got = reads_and_builds(&m, &part, &specs, oracle, 30.0, None);
        refused_at(&m, &part, &specs, got * 1.01);
    }
}

#[test]
fn no_specs_have_no_largest_measure() {
    let mut m = Model::new();
    let (bar, _, _) = bar(&mut m);
    assert!(largest_chamfer_distance_all(&m, &bar, &[], 10.0, T).is_err());
}
