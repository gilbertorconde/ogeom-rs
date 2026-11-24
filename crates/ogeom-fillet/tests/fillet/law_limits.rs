//! Largest-size readings for variable-radius laws (ADR 0001, amendment
//! "largest-size readings for laws"): one value of a law freed, or the
//! whole law scaled, each confirmed against the build. For each case the
//! reading's limit `L` builds valid at `L (1 - tol)` and is refused at
//! `L (1 + tol)`, `tol` the reading's own [`LIMIT_TOLERANCE`], and `L`
//! stands at the closed form the geometry gives.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{
    BlendLimit, BlendSizing, DistanceLawSize, FilletCorner, LIMIT_TOLERANCE, LawSize, LawStation,
    RadiusLaw, RoundProfile, blend_limit,
};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A request under a law, owned, and its free value.
enum Ask {
    Law(Vec<Shape>, RadiusLaw, LawSize),
    Linear(Shape, f64, f64, LawSize),
    Distances(Shape, Shape, RadiusLaw, RadiusLaw, DistanceLawSize),
}

/// `law`'s points, a linear law's two ends.
fn points(law: &RadiusLaw) -> Vec<(f64, f64)> {
    match law {
        RadiusLaw::Linear { start, end } => vec![(0.0, *start), (1.0, *end)],
        RadiusLaw::Points(p) => p.clone(),
    }
}

/// `law` with point `i` at `v` (both ends where `joined`), or scaled.
fn with(law: &RadiusLaw, free: Option<LawSize>, v: f64, joined: bool) -> RadiusLaw {
    let mut p = points(law);
    let n = p.len();
    match free {
        None => {}
        Some(LawSize::Scale) => p.iter_mut().for_each(|q| q.1 *= v),
        Some(LawSize::At(at)) => {
            let i = match at {
                LawStation::Start => 0,
                LawStation::End => n - 1,
                LawStation::Point(i) => i,
            };
            p[i].1 = v;
            if joined && (i == 0 || i == n - 1) {
                p[0].1 = v;
                p[n - 1].1 = v;
            }
        }
    }
    match law {
        RadiusLaw::Linear { .. } => RadiusLaw::Linear {
            start: p[0].1,
            end: p[1].1,
        },
        RadiusLaw::Points(_) => RadiusLaw::Points(p),
    }
}

impl Ask {
    fn sizing(&self) -> BlendSizing<'_> {
        match self {
            Self::Law(chain, law, size) => BlendSizing::RoundLaw {
                chain,
                law,
                size: *size,
            },
            Self::Linear(edge, start, end, size) => BlendSizing::RoundLinear {
                edge,
                start: *start,
                end: *end,
                size: *size,
            },
            Self::Distances(edge, face, on_face, on_other, size) => BlendSizing::DistanceLaws {
                edge,
                face,
                on_face,
                on_other,
                profile: RoundProfile::Elliptical,
                size: *size,
            },
        }
    }

    /// The public call at `v`, on a copy of the model; `joined` where the
    /// law's first and last points are one point of a closed chain.
    fn build(&self, m: &Model, s: &Shape, v: f64, joined: bool) -> Result<(Model, Shape), String> {
        let mut m = m.clone();
        let built = match self {
            Self::Law(chain, law, size) => ogeom_fillet::fillet_edges_variable(
                &mut m,
                s,
                chain,
                &with(law, Some(*size), v, joined),
                T,
            ),
            Self::Linear(edge, start, end, size) => {
                let RadiusLaw::Linear { start, end } = with(
                    &RadiusLaw::Linear {
                        start: *start,
                        end: *end,
                    },
                    Some(*size),
                    v,
                    false,
                ) else {
                    unreachable!()
                };
                ogeom_fillet::fillet_edge_variable(&mut m, s, edge, start, end, T)
            }
            Self::Distances(edge, face, on_face, on_other, size) => {
                let (a, b) = match size {
                    DistanceLawSize::OnFace(at) => (Some(LawSize::At(*at)), None),
                    DistanceLawSize::OnOther(at) => (None, Some(LawSize::At(*at))),
                    DistanceLawSize::Scale => (Some(LawSize::Scale), Some(LawSize::Scale)),
                };
                ogeom_fillet::fillet_edge_distances_variable(
                    &mut m,
                    s,
                    edge,
                    face,
                    &with(on_face, a, v, joined),
                    &with(on_other, b, v, joined),
                    RoundProfile::Elliptical,
                    T,
                )
            }
        };
        built.map(|b| (m, b.shape)).map_err(|e| e.to_string())
    }
}

/// The reading for `ask`, held to its promise: valid at `L (1 - tol)`,
/// refused at `L (1 + tol)`; the refusal's words returned.
fn confirmed(
    name: &str,
    m: &Model,
    s: &Shape,
    ask: &Ask,
    up_to: f64,
    joined: bool,
) -> (BlendLimit, String) {
    let found =
        blend_limit(m, s, &ask.sizing(), up_to, T).unwrap_or_else(|e| panic!("{name}: {e}"));
    let l = found
        .size
        .unwrap_or_else(|| panic!("{name}: no size builds: {found:?}"));
    let (built_m, built) = ask
        .build(m, s, l * (1.0 - LIMIT_TOLERANCE), joined)
        .unwrap_or_else(|e| panic!("{name}: refused below the limit {l}: {e}"));
    let diagnosis = ogeom_algo::check(&built_m, &built, T).unwrap();
    assert!(
        diagnosis.is_valid(),
        "{name}: invalid below the limit {l}: {:?}",
        diagnosis.problems
    );
    assert!(
        found.refused.is_some(),
        "{name}: the bound builds: {found:?}"
    );
    let refusal = ask
        .build(m, s, l * (1.0 + LIMIT_TOLERANCE), joined)
        .err()
        .unwrap_or_else(|| panic!("{name}: builds a tolerance past the limit {l}: {found:?}"));
    assert_eq!(found.builds_again, None, "{name}: {found:?}");
    (found, refusal)
}

/// `size` stands within the reading's tolerance below `exact`.
fn near(name: &str, size: f64, exact: f64) {
    assert!(
        size <= exact * (1.0 + LIMIT_TOLERANCE) && size >= exact * (1.0 - 2.0 * LIMIT_TOLERANCE),
        "{name}: limit {size} against {exact}"
    );
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

/// A prism of the (y, z) profile swept 30 along x: its cap at x = 0 is the
/// profile itself.
fn prism(profile: &[(f64, f64)]) -> (Model, Shape) {
    let mut m = Model::new();
    let points: Vec<Point> = profile
        .iter()
        .map(|&(y, z)| Point::new(0.0, y, z))
        .collect();
    let wire = ogeom_algo::make_polygon(&mut m, &points, true, T)
        .unwrap()
        .shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::X),
        (-1000.0, 1000.0),
        (-1000.0, 1000.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(&mut m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let swept = ogeom_algo::make_prism(&mut m, &face, Vector::new(30.0, 0.0, 0.0), T)
        .unwrap()
        .shape;
    let s = ogeom_algo::baked_shape(&mut m, &swept, T).unwrap().shape;
    (m, s)
}

/// A wedge whose cap at x = 0 deepens from 4 below its top edge at y = 20
/// to 6 at y = 0: the bottom slopes. The top edge on the cap (y from 0 to
/// 20 at z = 6) meets the cap square, so a ball of radius `r` touches the
/// cap `r` below the edge, and fits where `r` is at most the cap's depth
/// there, `6 - y / 10`.
fn wedge() -> (Model, Shape, Shape, Shape) {
    let (m, s) = prism(&[(0.0, 0.0), (20.0, 2.0), (20.0, 6.0), (0.0, 6.0)]);
    let edge = edge_near(&m, &s, [0.0, 10.0, 6.0]);
    let cap = face_near(&m, &s, [0.0, 10.0, 3.5]);
    (m, s, edge, cap)
}

/// The cap's depth below the top edge at the end the edge runs to.
fn depth_at_end(m: &Model, edge: &Shape) -> f64 {
    let (_, b) = ogeom_algo::edge_vertices(m, edge).unwrap().unwrap();
    let y = m.node(&b).unwrap().data().as_vertex().unwrap().point.y;
    6.0 - y / 10.0
}

/// A linear law on a straight edge between planes, its end radius free:
/// the limit is the cap's depth at that end, 4 or 6 as the edge runs,
/// the start held at 1. The ball there touches the cap's bottom corner;
/// one a tolerance larger leaves the cap.
#[test]
fn a_linear_law_s_end_radius_reaches_the_face_s_width_at_that_end() {
    let (m, s, edge, _) = wedge();
    let exact = depth_at_end(&m, &edge);
    let law = RadiusLaw::Linear {
        start: 1.0,
        end: 1.0,
    };
    let ask = Ask::Law(vec![edge.clone()], law, LawSize::At(LawStation::End));
    let (found, why) = confirmed("wedge law end", &m, &s, &ask, 20.0, false);
    near("wedge law end", found.size.unwrap(), exact);
    assert!(
        found.probes <= 10,
        "the seat brackets it: {found:?} ({why})"
    );

    let ask = Ask::Linear(edge, 1.0, 1.0, LawSize::At(LawStation::End));
    let (found, _) = confirmed("wedge linear end", &m, &s, &ask, 20.0, false);
    near("wedge linear end", found.size.unwrap(), exact);
}

/// A slab's top edge over a cap notched from below in the middle: the cap
/// is 5 deep but 2 deep under the edge's midpoint, the notch's apex at
/// (10, 3). A law of 1 at both ends and a free point halfway: the limit is
/// the narrow stretch's 2, where the ball halfway touches the apex; and
/// so is a distance law's along the cap. The
/// law peaks there with zero slope and the cap deepens either side, so no
/// other station comes nearer its face.
#[test]
fn an_interior_point_beside_a_narrow_face_reaches_that_face() {
    let (m, s) = prism(&[
        (0.0, 0.0),
        (8.0, 0.0),
        (10.0, 3.0),
        (12.0, 0.0),
        (20.0, 0.0),
        (20.0, 5.0),
        (0.0, 5.0),
    ]);
    let edge = edge_near(&m, &s, [0.0, 10.0, 5.0]);
    let law = RadiusLaw::Points(vec![(0.0, 1.0), (0.5, 1.0), (1.0, 1.0)]);
    let ask = Ask::Law(
        vec![edge.clone()],
        law.clone(),
        LawSize::At(LawStation::Point(1)),
    );
    let (found, _) = confirmed("notched interior point", &m, &s, &ask, 10.0, false);
    near("notched interior point", found.size.unwrap(), 2.0);

    // Two distance laws, the one along the cap free halfway, the one along
    // the top held at 1: the distance down the cap reaches the apex too.
    let cap = face_near(&m, &s, [0.0, 10.0, 2.5]);
    let one = RadiusLaw::Linear {
        start: 1.0,
        end: 1.0,
    };
    let ask = Ask::Distances(
        edge,
        cap,
        law,
        one,
        DistanceLawSize::OnFace(LawStation::Point(1)),
    );
    let (found, _) = confirmed("notched interior distance", &m, &s, &ask, 10.0, false);
    near("notched interior distance", found.size.unwrap(), 2.0);
}

/// A drum's top rim, a closed chain, under a law whose first and last
/// points are the rim's start, given twice: freeing the start frees both,
/// and the limit is the wall's height of 3, where the ball's contact
/// reaches the bottom rim. Read with only the first point moved, the law
/// would break at the start and no radius but the end's would build.
#[test]
fn a_closed_chain_s_law_wraps_round_its_start() {
    let mut m = Model::new();
    let s = ogeom_algo::make_cylinder(&mut m, at(0.0, 0.0, 0.0), 10.0, 3.0, T)
        .unwrap()
        .shape;
    let rim = edge_near(&m, &s, [0.0, -10.0, 3.0]);
    let law = RadiusLaw::Points(vec![(0.0, 1.0), (0.5, 1.0), (1.0, 1.0)]);
    let ask = Ask::Law(vec![rim], law, LawSize::At(LawStation::Start));
    let (found, _) = confirmed("drum rim wrapped start", &m, &s, &ask, 8.0, true);
    near("drum rim wrapped start", found.size.unwrap(), 3.0);
}

/// A law holding one radius scales to the constant round's limit over
/// that radius: a 2 mm plate's edge, whose round reaches the plate's
/// thickness, under a law of 0.5 scales to 4; and two distance laws of 3
/// and 1 scale as the constant pair does.
#[test]
fn a_constant_law_scales_to_the_constant_round_s_limit() {
    let mut m = Model::new();
    let s = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (20.0, 20.0, 2.0), T)
        .unwrap()
        .shape;
    let edge = edge_near(&m, &s, [10.0, 0.0, 2.0]);
    let top = face_near(&m, &s, [10.0, 10.0, 2.0]);

    let constant = blend_limit(
        &m,
        &s,
        &BlendSizing::Rounds {
            edges: std::slice::from_ref(&edge),
            corner: FilletCorner::Setback,
        },
        10.0,
        T,
    )
    .unwrap()
    .size
    .unwrap();
    let law = RadiusLaw::Linear {
        start: 0.5,
        end: 0.5,
    };
    let ask = Ask::Law(vec![edge.clone()], law, LawSize::Scale);
    let (found, _) = confirmed("plate law scale", &m, &s, &ask, 40.0, false);
    near("plate law scale", found.size.unwrap() * 0.5, constant);
    near("plate law scale against 2", constant, 2.0);

    let pair = blend_limit(
        &m,
        &s,
        &BlendSizing::RoundDistances {
            edge: &edge,
            face: &top,
            on_face: 3.0,
            on_other: 1.0,
            profile: RoundProfile::Elliptical,
        },
        30.0,
        T,
    )
    .unwrap()
    .size
    .unwrap();
    let ask = Ask::Distances(
        edge,
        top,
        RadiusLaw::Points(vec![(0.0, 3.0)]),
        RadiusLaw::Points(vec![(0.0, 1.0)]),
        DistanceLawSize::Scale,
    );
    let (found, _) = confirmed("plate distance laws scale", &m, &s, &ask, 30.0, false);
    near("plate distance laws scale", found.size.unwrap(), pair);
}

/// A point the law does not have is refused by name, before any build.
#[test]
fn a_point_the_law_lacks_is_refused() {
    let (m, s, edge, _) = wedge();
    let law = RadiusLaw::Points(vec![(0.0, 1.0), (1.0, 1.0)]);
    let err =
        ogeom_fillet::largest_law_radius(&m, &s, &[edge], &law, LawStation::Point(2), 10.0, T)
            .unwrap_err();
    assert!(err.to_string().contains("no point 2"), "{err}");
}
