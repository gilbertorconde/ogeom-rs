//! `BlendRequest` against the calls it routes to: for each row of its
//! dispatch, the request and the old call, each on its own copy of the same
//! part, give the same result (faces, area and volume to the bit) or the
//! same refusal text; and each combination no call builds is refused by
//! name (docs/design/1265-fillet-request.md).
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_fillet::{
    Blend, BlendLimit, BlendRefusal, BlendRequest, BlendSection, BlendSizing, Chamfer,
    ChamferCorner, DistanceLawSize, FilletCorner, LawSize, LawStation, Overflow, RadiusLaw,
    RoundProfile, Setbacks, Vary,
};
use ogeom_geom::{Continuity, Curve3d as _};
use ogeom_math::{Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The block 20 by 20 by 10 from the origin.
fn block() -> (Model, Shape) {
    let mut m = Model::new();
    let s = ogeom_algo::make_box(&mut m, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    (m, s)
}

/// The block 30 by 20 by 15 under `z = 0`, a pocket of radius 5 and depth
/// 8 sunk into its top about `(15, 7)`: a round of 4 on the front top edge
/// rolls over the pocket's rim, which hold edge keeps and which, with hold
/// edge and trim through off, is refused.
fn pocketed() -> (Model, Shape) {
    let mut m = Model::new();
    let at = |x: f64, y: f64, z: f64| {
        Frame::new(
            Point::new(x, y, z),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap()
    };
    let block = ogeom_algo::make_box(&mut m, at(0.0, 0.0, -15.0), (30.0, 20.0, 15.0), T)
        .unwrap()
        .shape;
    let pocket = ogeom_algo::make_cylinder(&mut m, at(15.0, 7.0, -8.0), 5.0, 9.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(&mut m, &block, &pocket, T).unwrap().shape;
    (m, s)
}

/// The pocketed block's front top edge, alone in a one-edge list.
fn foot(m: &Model, s: &Shape) -> Vec<Shape> {
    let found = edges_where(m, s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    assert_eq!(found.len(), 1, "the front top edge");
    found
}

fn edges_where(m: &Model, s: &Shape, on: impl Fn(Point) -> bool) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            [0.0, 0.3, 0.7, 1.0]
                .iter()
                .all(|t| on(c.point_at(a + (b - a) * t, T).unwrap()))
        })
        .collect()
}

/// The block's edge through `p`.
fn edge(m: &Model, s: &Shape, p: [f64; 3]) -> Shape {
    let q = Point::new(p[0], p[1], p[2]);
    let through: Vec<Shape> = explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            (0..=16).any(|i| {
                c.point_at(a + (b - a) * f64::from(i) / 16.0, T)
                    .unwrap()
                    .distance(q)
                    < 1e-9
            })
        })
        .collect();
    assert_eq!(through.len(), 1, "one edge through {p:?}");
    through[0].clone()
}

/// The three edges leaving the block's top corner `(0, 0, 10)`: along `x`,
/// along `y`, and down `z`.
fn corner(m: &Model, s: &Shape) -> Vec<Shape> {
    [[10.0, 0.0, 10.0], [0.0, 10.0, 10.0], [0.0, 0.0, 5.0]]
        .iter()
        .map(|p| edge(m, s, *p))
        .collect()
}

fn corner_vertex(m: &Model, s: &Shape) -> Shape {
    explore_unique(m, s, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .find(|v| {
            m.node(v)
                .and_then(|n| n.data().as_vertex().map(|d| d.point))
                .is_some_and(|x| x.distance(Point::new(0.0, 0.0, 10.0)) < 1e-9)
        })
        .expect("the corner vertex")
}

/// The front face `y = 0`, one of the upright edge's two faces.
fn front(m: &Model, s: &Shape, e: &Shape) -> Shape {
    ogeom_topo::ancestors_of(m, s, e, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let probe = edges_where(m, f, |p| p.y.abs() < 1e-9);
            probe.len() == 4
        })
        .expect("the front face")
}

fn faces(m: &Model, s: &Shape) -> usize {
    explore_unique(m, s, ShapeType::Face).unwrap().len()
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

fn area(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::surface_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// What a result is, as the comparison reads it: faces, then area and
/// volume by their bits; or the refusal's text.
fn read(m: &Model, built: &OgeomResult<Built>) -> Result<(usize, u64, u64), String> {
    match built {
        Ok(b) => Ok((
            faces(m, &b.shape),
            area(m, &b.shape).to_bits(),
            volume(m, &b.shape).to_bits(),
        )),
        Err(e) => Err(e.to_string()),
    }
}

type Call<'a, R> = &'a dyn Fn(&mut Model, &Shape) -> R;

/// `old` and `new` on their own copies of `part`, read alike; `built`
/// whether the old call is expected to build.
fn same(
    part: fn() -> (Model, Shape),
    built: bool,
    old: Call<'_, OgeomResult<Built>>,
    new: Call<'_, OgeomResult<Built>>,
) {
    let (mut a, s) = part();
    let x = old(&mut a, &s);
    let (mut b, t) = part();
    let y = new(&mut b, &t);
    let (x, y) = (read(&a, &x), read(&b, &y));
    assert_eq!(x.is_ok(), built, "the old call: {x:?}");
    assert_eq!(x, y);
}

/// As [`same`] for the checked reading: the refusal's edges and text.
fn same_checked(
    part: fn() -> (Model, Shape),
    built: bool,
    old: Call<'_, Result<Built, BlendRefusal>>,
    new: Call<'_, Result<Built, BlendRefusal>>,
) {
    let (mut a, s) = part();
    let x = old(&mut a, &s);
    let (mut b, t) = part();
    let y = new(&mut b, &t);
    let read_checked = |m: &Model, r: Result<Built, BlendRefusal>| match r {
        Ok(b) => Ok(read(m, &Ok(b)).unwrap()),
        Err(refusal) => Err((refusal.edges, refusal.error.to_string())),
    };
    let (x, y) = (read_checked(&a, x), read_checked(&b, y));
    assert_eq!(x.is_ok(), built, "the old call: {x:?}");
    assert_eq!(x, y);
}

/// As [`same`] for the partial reading: the result and every refusal.
fn same_partial(
    part: fn() -> (Model, Shape),
    old: Call<'_, OgeomResult<(Built, Vec<BlendRefusal>)>>,
    new: Call<'_, OgeomResult<(Built, Vec<BlendRefusal>)>>,
) -> usize {
    let (mut a, s) = part();
    let x = old(&mut a, &s).unwrap();
    let (mut b, t) = part();
    let y = new(&mut b, &t).unwrap();
    let refusals = |r: &[BlendRefusal]| -> Vec<(Vec<usize>, String)> {
        r.iter()
            .map(|r| (r.edges.clone(), r.error.to_string()))
            .collect()
    };
    assert_eq!(read(&a, &Ok(x.0)), read(&b, &Ok(y.0)));
    assert_eq!(refusals(&x.1), refusals(&y.1));
    x.1.len()
}

fn with(f: impl FnOnce(&mut Overflow)) -> Overflow {
    let mut o = Overflow::default();
    f(&mut o);
    o
}

fn no_hold_edge() -> Overflow {
    with(|o| {
        o.hold_edge = false;
        o.trim_through = false;
    })
}

fn hump() -> RadiusLaw {
    RadiusLaw::Points(vec![(0.0, 1.0), (0.5, 2.0), (1.0, 1.0)])
}

fn linear(start: f64, end: f64) -> RadiusLaw {
    RadiusLaw::Linear { start, end }
}

fn even(setback: f64) -> Setbacks {
    Setbacks {
        per_edge: Vec::new(),
        default: setback,
        continuity: Continuity::G1,
    }
}

// Rounds.

#[test]
fn rounds_of_one_radius_are_fillet_edges_overflow() {
    same(
        block,
        true,
        &|m, s| ogeom_fillet::fillet_edges(m, s, &corner(m, s), 2.0, T),
        &|m, s| BlendRequest::new().round(corner(m, s), 2.0).build(m, s, T),
    );
    for (overflow, built) in [(Overflow::default(), true), (no_hold_edge(), false)] {
        same(
            pocketed,
            built,
            &|m, s| ogeom_fillet::fillet_edges_overflow(m, s, &foot(m, s), 4.0, overflow, T),
            &|m, s| {
                BlendRequest::new()
                    .round(foot(m, s), 4.0)
                    .overflow(overflow)
                    .build(m, s, T)
            },
        );
    }
}

#[test]
fn rolling_rounds_are_fillet_edges_corner_overflow() {
    for overflow in [Overflow::default(), Overflow::none()] {
        same(
            block,
            true,
            &|m, s| {
                ogeom_fillet::fillet_edges_corner_overflow(
                    m,
                    s,
                    &corner(m, s),
                    2.0,
                    FilletCorner::Rolling,
                    overflow,
                    T,
                )
            },
            &|m, s| {
                BlendRequest::new()
                    .round(corner(m, s), 2.0)
                    .round_corner(FilletCorner::Rolling)
                    .overflow(overflow)
                    .build(m, s, T)
            },
        );
    }
}

#[test]
fn rounds_with_setbacks_are_fillet_edges_cornered() {
    same(
        block,
        true,
        &|m, s| {
            let v = corner_vertex(m, s);
            ogeom_fillet::fillet_edges_cornered(m, s, &corner(m, s), 2.0, &[(v, even(4.0))], T)
        },
        &|m, s| {
            BlendRequest::new()
                .round(corner(m, s), 2.0)
                .setback(corner_vertex(m, s), even(4.0))
                .build(m, s, T)
        },
    );
}

#[test]
fn rounds_of_several_radii_are_fillet_edges_with_or_blend_edges_overflow() {
    let radii = [1.0, 2.0, 1.5];
    let pairs = |m: &Model, s: &Shape| -> Vec<(Shape, f64)> {
        corner(m, s).into_iter().zip(radii).collect()
    };
    let request = |m: &Model, s: &Shape| {
        pairs(m, s)
            .into_iter()
            .fold(BlendRequest::new(), |r, (e, radius)| {
                r.edge(e, BlendSection::Round(radius))
            })
    };
    same(
        block,
        true,
        &|m, s| ogeom_fillet::fillet_edges_with(m, s, &pairs(m, s), T),
        &|m, s| request(m, s).build(m, s, T),
    );
    let overflow = Overflow::none();
    same(
        block,
        true,
        &|m, s| {
            let blends: Vec<(Shape, Blend)> = pairs(m, s)
                .into_iter()
                .map(|(e, r)| (e, Blend::Round(r)))
                .collect();
            ogeom_fillet::blend_edges_overflow(m, s, &blends, ChamferCorner::Facet, overflow, T)
        },
        &|m, s| request(m, s).overflow(overflow).build(m, s, T),
    );
}

/// A request's empty calls, a radius that rounds nothing and an edge asked
/// twice take the old call's refusal word for word.
#[test]
fn refusals_are_the_old_calls_own() {
    same(
        block,
        false,
        &|m, s| ogeom_fillet::fillet_edges(m, s, &[], 2.0, T),
        &|m, s| BlendRequest::new().build(m, s, T),
    );
    same(
        block,
        false,
        &|m, s| ogeom_fillet::fillet_edges(m, s, &[], 2.0, T),
        &|m, s| BlendRequest::new().round([], 2.0).build(m, s, T),
    );
    same(
        block,
        false,
        &|m, s| ogeom_fillet::chamfer_edges(m, s, &[], 1.0, T),
        &|m, s| BlendRequest::new().chamfer([], 1.0).build(m, s, T),
    );
    same(
        block,
        false,
        &|m, s| ogeom_fillet::fillet_edges_variable(m, s, &[], &hump(), T),
        &|m, s| {
            BlendRequest::new()
                .chain([], BlendSection::RoundLaw(hump()))
                .build(m, s, T)
        },
    );
    same(
        block,
        false,
        &|m, s| ogeom_fillet::fillet_edges(m, s, &corner(m, s), 0.0, T),
        &|m, s| BlendRequest::new().round(corner(m, s), 0.0).build(m, s, T),
    );
    same(
        block,
        false,
        &|m, s| {
            let e = corner(m, s)[0].clone();
            ogeom_fillet::chamfer_edges(m, s, &[e.clone(), e], 1.0, T)
        },
        &|m, s| {
            let e = corner(m, s)[0].clone();
            BlendRequest::new()
                .chamfer([e.clone(), e], 1.0)
                .build(m, s, T)
        },
    );
}

// Laws and distances.

#[test]
fn a_round_law_is_fillet_edges_variable_overflow() {
    for overflow in [Overflow::default(), Overflow::none()] {
        same(
            block,
            true,
            &|m, s| {
                let e = corner(m, s)[2].clone();
                ogeom_fillet::fillet_edges_variable_overflow(m, s, &[e], &hump(), overflow, T)
            },
            &|m, s| {
                let e = corner(m, s)[2].clone();
                BlendRequest::new()
                    .chain([e], BlendSection::RoundLaw(hump()))
                    .overflow(overflow)
                    .build(m, s, T)
            },
        );
    }
}

#[test]
fn a_round_of_two_distances_is_fillet_edge_distances_overflow() {
    for profile in [RoundProfile::Elliptical, RoundProfile::Conic(0.6)] {
        same(
            block,
            true,
            &|m, s| {
                let e = corner(m, s)[2].clone();
                let f = front(m, s, &e);
                ogeom_fillet::fillet_edge_distances(m, s, &e, &f, 1.0, 2.0, profile, T)
            },
            &|m, s| {
                let e = corner(m, s)[2].clone();
                let face = front(m, s, &e);
                BlendRequest::new()
                    .edge(
                        e,
                        BlendSection::RoundDistances {
                            face,
                            on_face: 1.0,
                            on_other: 2.0,
                            profile,
                        },
                    )
                    .build(m, s, T)
            },
        );
    }
}

#[test]
fn distance_laws_are_fillet_edges_distances_variable() {
    same(
        block,
        true,
        &|m, s| {
            let e = corner(m, s)[2].clone();
            let f = front(m, s, &e);
            ogeom_fillet::fillet_edges_distances_variable(
                m,
                s,
                &[e],
                &f,
                &linear(1.0, 2.0),
                &linear(1.5, 1.0),
                RoundProfile::Elliptical,
                T,
            )
        },
        &|m, s| {
            let e = corner(m, s)[2].clone();
            let face = front(m, s, &e);
            BlendRequest::new()
                .chain(
                    [e],
                    BlendSection::DistanceLaws {
                        face,
                        on_face: linear(1.0, 2.0),
                        on_other: linear(1.5, 1.0),
                        profile: RoundProfile::Elliptical,
                    },
                )
                .build(m, s, T)
        },
    );
}

// Chamfers and mixed requests.

fn specs(m: &Model, s: &Shape) -> Vec<(Shape, Chamfer)> {
    let edges = corner(m, s);
    let face = front(m, s, &edges[2]);
    vec![
        (edges[0].clone(), Chamfer::Symmetric(1.0)),
        (edges[1].clone(), Chamfer::Symmetric(1.5)),
        (
            edges[2].clone(),
            Chamfer::Distances {
                face,
                on_face: 1.0,
                on_other: 2.0,
            },
        ),
    ]
}

fn bevels(m: &Model, s: &Shape) -> BlendRequest {
    specs(m, s)
        .into_iter()
        .fold(BlendRequest::new(), |r, (e, c)| r.edge(e, c))
}

#[test]
fn bevels_at_the_defaults_are_chamfer_edges_with() {
    same(
        block,
        true,
        &|m, s| ogeom_fillet::chamfer_edges_with(m, s, &specs(m, s), T),
        &|m, s| bevels(m, s).build(m, s, T),
    );
    same(
        block,
        true,
        &|m, s| ogeom_fillet::chamfer_edges(m, s, &corner(m, s), 1.0, T),
        &|m, s| {
            BlendRequest::new()
                .chamfer(corner(m, s), 1.0)
                .build(m, s, T)
        },
    );
}

#[test]
fn bevels_under_other_options_are_chamfer_edges_overflow() {
    for (corner_kind, overflow) in [
        (ChamferCorner::Point, Overflow::default()),
        (ChamferCorner::Point, Overflow::none()),
        (ChamferCorner::Facet, Overflow::none()),
    ] {
        same(
            block,
            true,
            &|m, s| {
                ogeom_fillet::chamfer_edges_overflow(m, s, &specs(m, s), corner_kind, overflow, T)
            },
            &|m, s| {
                bevels(m, s)
                    .chamfer_corner(corner_kind)
                    .overflow(overflow)
                    .build(m, s, T)
            },
        );
    }
    same(
        block,
        true,
        &|m, s| ogeom_fillet::chamfer_edges_cornered(m, s, &specs(m, s), ChamferCorner::Point, T),
        &|m, s| {
            bevels(m, s)
                .chamfer_corner(ChamferCorner::Point)
                .build(m, s, T)
        },
    );
}

fn mixed(m: &Model, s: &Shape) -> Vec<(Shape, Blend)> {
    let edges = corner(m, s);
    let face = front(m, s, &edges[2]);
    vec![
        (edges[0].clone(), Blend::Round(2.0)),
        (edges[1].clone(), Blend::Chamfer(Chamfer::Symmetric(1.5))),
        (
            edges[2].clone(),
            Blend::RoundDistances {
                face,
                on_face: 1.0,
                on_other: 2.0,
                profile: RoundProfile::Elliptical,
            },
        ),
    ]
}

#[test]
fn blends_of_several_kinds_are_blend_edges_overflow() {
    let request = |m: &Model, s: &Shape| {
        mixed(m, s)
            .into_iter()
            .fold(BlendRequest::new(), |r, (e, b)| r.edge(e, b))
    };
    for corner_kind in [ChamferCorner::Facet, ChamferCorner::Point] {
        same(
            block,
            true,
            &|m, s| ogeom_fillet::blend_edges(m, s, &mixed(m, s), corner_kind, T),
            &|m, s| request(m, s).chamfer_corner(corner_kind).build(m, s, T),
        );
    }
    let two: fn(&Model, &Shape) -> Vec<(Shape, Blend)> = |m, s| mixed(m, s)[..2].to_vec();
    same(
        block,
        true,
        &|m, s| {
            ogeom_fillet::blend_edges_overflow(
                m,
                s,
                &two(m, s),
                ChamferCorner::Facet,
                Overflow::none(),
                T,
            )
        },
        &|m, s| {
            two(m, s)
                .into_iter()
                .fold(BlendRequest::new(), |r, (e, b)| r.edge(e, b))
                .overflow(Overflow::none())
                .build(m, s, T)
        },
    );
}

// Checked and partial readings.

#[test]
fn checked_rounds_and_bevels_are_the_in_order_readings() {
    for (radius, built) in [(2.0, true), (12.0, false)] {
        for corner_kind in [FilletCorner::Setback, FilletCorner::Rolling] {
            same_checked(
                block,
                built,
                &|m, s| {
                    ogeom_fillet::fillet_edges_corner_overflow_checked(
                        m,
                        s,
                        &corner(m, s),
                        radius,
                        corner_kind,
                        Overflow::default(),
                        T,
                    )
                },
                &|m, s| {
                    BlendRequest::new()
                        .round(corner(m, s), radius)
                        .round_corner(corner_kind)
                        .build_checked(m, s, T)
                },
            );
        }
    }
    for overflow in [Overflow::default(), no_hold_edge()] {
        same_checked(
            pocketed,
            overflow == Overflow::default(),
            &|m, s| {
                ogeom_fillet::fillet_edges_corner_overflow_checked(
                    m,
                    s,
                    &foot(m, s),
                    4.0,
                    FilletCorner::Setback,
                    overflow,
                    T,
                )
            },
            &|m, s| {
                BlendRequest::new()
                    .round(foot(m, s), 4.0)
                    .overflow(overflow)
                    .build_checked(m, s, T)
            },
        );
    }
    for (d, built) in [(1.0, true), (12.0, false)] {
        same_checked(
            block,
            built,
            &|m, s| {
                let specs: Vec<(Shape, Chamfer)> = corner(m, s)
                    .into_iter()
                    .map(|e| (e, Chamfer::Symmetric(d)))
                    .collect();
                ogeom_fillet::chamfer_edges_overflow_checked(
                    m,
                    s,
                    &specs,
                    ChamferCorner::Point,
                    Overflow::none(),
                    T,
                )
            },
            &|m, s| {
                BlendRequest::new()
                    .chamfer(corner(m, s), d)
                    .chamfer_corner(ChamferCorner::Point)
                    .overflow(Overflow::none())
                    .build_checked(m, s, T)
            },
        );
    }
}

#[test]
fn partial_rounds_and_bevels_are_the_partial_readings() {
    // Rounds of 6 on the corner's edges and one more on the far top edge:
    // the upright edge is 10 tall, so the corner refuses.
    let edges = |m: &Model, s: &Shape| {
        let mut e = corner(m, s);
        e.push(edge(m, s, [20.0, 10.0, 10.0]));
        e
    };
    for radius in [2.0, 6.0] {
        same_partial(
            block,
            &|m, s| ogeom_fillet::fillet_edges_partial(m, s, &edges(m, s), radius, T),
            &|m, s| {
                BlendRequest::new()
                    .round(edges(m, s), radius)
                    .build_partial(m, s, T)
            },
        );
    }
    same_partial(
        block,
        &|m, s| ogeom_fillet::chamfer_edges_partial(m, s, &specs(m, s), T),
        &|m, s| bevels(m, s).build_partial(m, s, T),
    );
}

// Largest-size readings.

fn same_limit(
    part: fn() -> (Model, Shape),
    old: &dyn Fn(&Model, &Shape) -> OgeomResult<BlendLimit>,
    new: &dyn Fn(&Model, &Shape) -> OgeomResult<BlendLimit>,
) -> f64 {
    let (a, s) = part();
    let (b, t) = part();
    let (x, y) = (old(&a, &s).unwrap(), new(&b, &t).unwrap());
    assert_eq!(format!("{x:?}"), format!("{y:?}"));
    x.size.expect("a size builds")
}

#[test]
fn round_readings_are_blend_limit_over_their_sizing() {
    let up_to = 10.0;
    let r = same_limit(
        block,
        &|m, s| {
            let edges = corner(m, s);
            let sizing = BlendSizing::Rounds {
                edges: &edges,
                corner: FilletCorner::Setback,
            };
            ogeom_fillet::blend_limit(m, s, &sizing, up_to, T)
        },
        &|m, s| {
            BlendRequest::new()
                .round(corner(m, s), 1.0)
                .limit(m, s, Vary::Scale, up_to, T)
        },
    );
    let (m, s) = block();
    let old = ogeom_fillet::largest_fillet_radius_all(&m, &s, &corner(&m, &s), up_to, T).unwrap();
    let new = BlendRequest::new()
        .round(corner(&m, &s), 1.0)
        .largest(&m, &s, Vary::Scale, up_to, T)
        .unwrap();
    assert_eq!(old.map(f64::to_bits), new.map(f64::to_bits));
    assert_eq!(old, Some(r));
    // A request at radius 2 reads the factor: the radius halved.
    let half = BlendRequest::new()
        .round(corner(&m, &s), 2.0)
        .largest(&m, &s, Vary::Scale, up_to / 2.0, T)
        .unwrap();
    assert_eq!(half, Some(r / 2.0));
    let old = ogeom_fillet::blend_seat(
        &m,
        &s,
        &BlendSizing::Rounds {
            edges: &corner(&m, &s),
            corner: FilletCorner::Setback,
        },
        up_to,
        T,
    )
    .unwrap();
    let new = BlendRequest::new()
        .round(corner(&m, &s), 1.0)
        .seat(&m, &s, Vary::Scale, up_to, T)
        .unwrap();
    assert_eq!(old.map(f64::to_bits), new.map(f64::to_bits));

    let overflow = Overflow::none();
    same_limit(
        block,
        &|m, s| {
            let edges = corner(m, s);
            let sizing = BlendSizing::RoundsOverflow {
                edges: &edges,
                overflow,
            };
            ogeom_fillet::blend_limit(m, s, &sizing, up_to, T)
        },
        &|m, s| {
            BlendRequest::new()
                .round(corner(m, s), 1.0)
                .overflow(overflow)
                .limit(m, s, Vary::Scale, up_to, T)
        },
    );
    same_limit(
        block,
        &|m, s| {
            let edges = corner(m, s);
            let sizing = BlendSizing::RoundsCornerOverflow {
                edges: &edges,
                corner: FilletCorner::Rolling,
                overflow,
            };
            ogeom_fillet::blend_limit(m, s, &sizing, up_to, T)
        },
        &|m, s| {
            BlendRequest::new()
                .round(corner(m, s), 1.0)
                .round_corner(FilletCorner::Rolling)
                .overflow(overflow)
                .limit(m, s, Vary::Scale, up_to, T)
        },
    );
}

#[test]
fn distance_and_bevel_readings_are_blend_limit_over_their_sizing() {
    let up_to = 8.0;
    let (m, s) = block();
    let e = corner(&m, &s)[2].clone();
    let face = front(&m, &s, &e);
    let old =
        ogeom_fillet::largest_fillet_distances(&m, &s, &e, &face, 1.0, 2.0, up_to, T).unwrap();
    let new = BlendRequest::new()
        .edge(
            e.clone(),
            BlendSection::RoundDistances {
                face: face.clone(),
                on_face: 1.0,
                on_other: 2.0,
                profile: RoundProfile::Elliptical,
            },
        )
        .largest(&m, &s, Vary::Scale, up_to, T)
        .unwrap();
    assert!(old.is_some());
    assert_eq!(old.map(f64::to_bits), new.map(f64::to_bits));

    for spec in [
        Chamfer::Symmetric(1.0),
        Chamfer::Angle {
            face: face.clone(),
            distance: 1.0,
            angle: 0.5,
        },
        Chamfer::Distances {
            face: face.clone(),
            on_face: 1.0,
            on_other: 2.0,
        },
    ] {
        let old = ogeom_fillet::largest_chamfer_distance(&m, &s, &e, &spec, up_to, T).unwrap();
        let new = BlendRequest::new()
            .edge(e.clone(), spec.clone())
            .largest(&m, &s, Vary::Scale, up_to, T)
            .unwrap();
        assert!(old.is_some());
        assert_eq!(old.map(f64::to_bits), new.map(f64::to_bits), "{spec:?}");
    }
    let all: Vec<(Shape, Chamfer)> = corner(&m, &s)
        .into_iter()
        .map(|e| (e, Chamfer::Symmetric(1.0)))
        .collect();
    let old = ogeom_fillet::largest_chamfer_distance_all(&m, &s, &all, up_to, T).unwrap();
    let new = BlendRequest::new()
        .chamfer(corner(&m, &s), 1.0)
        .largest(&m, &s, Vary::Scale, up_to, T)
        .unwrap();
    assert!(old.is_some());
    assert_eq!(old.map(f64::to_bits), new.map(f64::to_bits));
    let old = ogeom_fillet::largest_chamfer_distance_all_overflow(
        &m,
        &s,
        &all,
        up_to,
        ChamferCorner::Point,
        Overflow::none(),
        T,
    )
    .unwrap();
    let new = BlendRequest::new()
        .chamfer(corner(&m, &s), 1.0)
        .chamfer_corner(ChamferCorner::Point)
        .overflow(Overflow::none())
        .largest(&m, &s, Vary::Scale, up_to, T)
        .unwrap();
    assert_eq!(old.map(f64::to_bits), new.map(f64::to_bits));
}

#[test]
fn law_readings_are_blend_limit_over_their_sizing() {
    let up_to = 8.0;
    let (m, s) = block();
    let e = corner(&m, &s)[2].clone();
    let face = front(&m, &s, &e);
    let law = || BlendRequest::new().chain([e.clone()], BlendSection::RoundLaw(hump()));
    let old = ogeom_fillet::largest_law_radius(
        &m,
        &s,
        std::slice::from_ref(&e),
        &hump(),
        LawStation::Point(1),
        up_to,
        T,
    )
    .unwrap();
    let new = law()
        .largest(&m, &s, Vary::Law(LawStation::Point(1)), up_to, T)
        .unwrap();
    assert!(old.is_some());
    assert_eq!(old.map(f64::to_bits), new.map(f64::to_bits));
    let old = ogeom_fillet::largest_law_scale(&m, &s, std::slice::from_ref(&e), &hump(), up_to, T)
        .unwrap();
    let new = law().largest(&m, &s, Vary::Scale, up_to, T).unwrap();
    assert!(old.is_some());
    assert_eq!(old.map(f64::to_bits), new.map(f64::to_bits));
    let law_overflow = hump();
    let sizing = BlendSizing::RoundLawOverflow {
        chain: std::slice::from_ref(&e),
        law: &law_overflow,
        size: LawSize::Scale,
        overflow: Overflow::none(),
    };
    let old = ogeom_fillet::blend_limit(&m, &s, &sizing, up_to, T).unwrap();
    let new = law()
        .overflow(Overflow::none())
        .limit(&m, &s, Vary::Scale, up_to, T)
        .unwrap();
    assert_eq!(format!("{old:?}"), format!("{new:?}"));

    for (size, vary) in [
        (
            DistanceLawSize::OnFace(LawStation::Start),
            Vary::OnFace(LawStation::Start),
        ),
        (
            DistanceLawSize::OnOther(LawStation::End),
            Vary::OnOther(LawStation::End),
        ),
        (DistanceLawSize::Scale, Vary::Scale),
    ] {
        let old = ogeom_fillet::largest_law_distances(
            &m,
            &s,
            &e,
            &face,
            &linear(1.0, 2.0),
            &linear(1.5, 1.0),
            RoundProfile::Elliptical,
            size,
            up_to,
            T,
        )
        .unwrap();
        let new = BlendRequest::new()
            .edge(
                e.clone(),
                BlendSection::DistanceLaws {
                    face: face.clone(),
                    on_face: linear(1.0, 2.0),
                    on_other: linear(1.5, 1.0),
                    profile: RoundProfile::Elliptical,
                },
            )
            .largest(&m, &s, vary, up_to, T)
            .unwrap();
        assert!(old.is_some());
        assert_eq!(old.map(f64::to_bits), new.map(f64::to_bits), "{vary:?}");
    }
}

// Combinations no call builds, each refused by name.

fn refusal(r: OgeomResult<impl core::fmt::Debug>) -> String {
    r.unwrap_err().to_string()
}

#[test]
fn combinations_no_call_builds_are_refused_by_name() {
    let (m, s) = block();
    let edges = corner(&m, &s);
    let v = corner_vertex(&m, &s);
    let e = edges[2].clone();
    let face = front(&m, &s, &e);
    let law = BlendSection::RoundLaw(hump());
    let distances = BlendSection::RoundDistances {
        face: face.clone(),
        on_face: 1.0,
        on_other: 2.0,
        profile: RoundProfile::Elliptical,
    };
    let distance_laws = BlendSection::DistanceLaws {
        face,
        on_face: linear(1.0, 2.0),
        on_other: linear(1.5, 1.0),
        profile: RoundProfile::Elliptical,
    };
    let several = BlendRequest::new()
        .edge(edges[0].clone(), BlendSection::Round(1.0))
        .edge(edges[1].clone(), BlendSection::Round(2.0));
    let mix = BlendRequest::new()
        .round([edges[0].clone()], 1.0)
        .chamfer([edges[1].clone()], 1.0);
    let rolling = |r: BlendRequest| r.round_corner(FilletCorner::Rolling);
    let build = |r: &BlendRequest| {
        let mut m = m.clone();
        refusal(r.build(&mut m, &s, T))
    };
    let checked = |r: &BlendRequest| {
        let mut m = m.clone();
        let refused = r.build_checked(&mut m, &s, T).unwrap_err();
        assert!(refused.edges.is_empty(), "{refused}");
        refused.error.to_string()
    };
    let partial = |r: &BlendRequest| {
        let mut m = m.clone();
        refusal(r.build_partial(&mut m, &s, T))
    };
    let reading = |r: &BlendRequest, vary: Vary| refusal(r.limit(&m, &s, vary, 8.0, T));
    let refused = |text: &str| format!("construction failed: {text}");

    let cases: Vec<(String, &str)> = vec![
        (
            build(
                &BlendRequest::new()
                    .round([edges[0].clone()], 1.0)
                    .chain([e.clone()], law.clone()),
            ),
            "a law with other blends in one request is not built",
        ),
        (
            build(
                &BlendRequest::new()
                    .chain([e.clone()], law.clone())
                    .chain([edges[0].clone()], law.clone()),
            ),
            "a law with other blends in one request is not built",
        ),
        (
            build(&rolling(
                BlendRequest::new().chain([e.clone()], law.clone()),
            )),
            "rolling corners on a round law are not built",
        ),
        (
            build(&rolling(several.clone())),
            "rolling corners on rounds of several radii are not built",
        ),
        (
            build(&rolling(
                BlendRequest::new().edge(e.clone(), distances.clone()),
            )),
            "rolling corners on a round of two distances are not built",
        ),
        (
            build(&rolling(mix.clone())),
            "rolling corners on blends of several kinds are not built",
        ),
        (
            build(&rolling(
                BlendRequest::new().chain([e.clone()], distance_laws.clone()),
            )),
            "rolling corners on distance laws are not built",
        ),
        (
            build(
                &BlendRequest::new()
                    .chamfer(edges.clone(), 1.0)
                    .setback(v.clone(), even(4.0)),
            ),
            "setback corners at vertices on bevels are not built",
        ),
        (
            build(&several.clone().setback(v.clone(), even(4.0))),
            "setback corners at vertices on rounds of several radii are not built",
        ),
        (
            build(
                &rolling(BlendRequest::new().round(edges.clone(), 2.0))
                    .setback(v.clone(), even(4.0)),
            ),
            "setback corners at vertices with rolling corners are not built",
        ),
        (
            build(
                &BlendRequest::new()
                    .round(edges.clone(), 2.0)
                    .setback(v.clone(), even(4.0))
                    .overflow(Overflow::none()),
            ),
            "setback corners at vertices under overflow options are not built",
        ),
        (
            build(
                &BlendRequest::new()
                    .chain([e.clone()], distance_laws.clone())
                    .overflow(Overflow::none()),
            ),
            "distance laws under overflow options are not built",
        ),
        (
            build(&BlendRequest::new().chain([], distances.clone())),
            "a round of two distances on no edge rounds nothing",
        ),
        (
            checked(&several),
            "a checked reading of rounds of several radii is not built",
        ),
        (
            checked(&mix),
            "a checked reading of blends of several kinds is not built",
        ),
        (
            checked(&BlendRequest::new().chain([e.clone()], law.clone())),
            "a checked reading of a round law is not built",
        ),
        (
            checked(
                &BlendRequest::new()
                    .round(edges.clone(), 2.0)
                    .setback(v.clone(), even(4.0)),
            ),
            "a checked reading with setback corners at vertices is not built",
        ),
        (
            partial(&several),
            "a partial reading of rounds of several radii is not built",
        ),
        (
            partial(
                &BlendRequest::new()
                    .round(edges.clone(), 2.0)
                    .overflow(Overflow::none()),
            ),
            "a partial reading under overflow options is not built",
        ),
        (
            partial(&rolling(BlendRequest::new().round(edges.clone(), 2.0))),
            "a partial reading with rolling corners is not built",
        ),
        (
            partial(
                &BlendRequest::new()
                    .chamfer(edges.clone(), 1.0)
                    .chamfer_corner(ChamferCorner::Point),
            ),
            "a partial reading with point corners is not built",
        ),
        (
            partial(
                &BlendRequest::new()
                    .round(edges.clone(), 2.0)
                    .setback(v.clone(), even(4.0)),
            ),
            "a partial reading with setback corners at vertices is not built",
        ),
        (
            reading(&several, Vary::Scale),
            "a reading of rounds of several radii is not read",
        ),
        (
            reading(&mix, Vary::Scale),
            "a reading of blends of several kinds is not read",
        ),
        (
            reading(
                &BlendRequest::new()
                    .round(edges.clone(), 1.0)
                    .setback(v.clone(), even(4.0)),
                Vary::Scale,
            ),
            "a reading of rounds with setback corners at vertices is not read",
        ),
        (
            reading(
                &BlendRequest::new().round(edges.clone(), 1.0),
                Vary::Law(LawStation::Start),
            ),
            "a law reading of rounds, a request with no law, is not read",
        ),
        (
            reading(
                &BlendRequest::new().chain([e.clone()], law.clone()),
                Vary::OnFace(LawStation::Start),
            ),
            "a distance-law reading of a round law is not read",
        ),
        (
            reading(
                &BlendRequest::new().chain([e.clone()], distance_laws.clone()),
                Vary::Law(LawStation::Start),
            ),
            "a round-law reading of distance laws is not read",
        ),
        (
            reading(
                &BlendRequest::new().chain(edges[..2].to_vec(), distance_laws.clone()),
                Vary::Scale,
            ),
            "a reading of distance laws over a chain of 2 edges is not read",
        ),
        (
            reading(
                &BlendRequest::new()
                    .chain([e.clone()], distance_laws.clone())
                    .overflow(Overflow::none()),
                Vary::Scale,
            ),
            "a reading of distance laws under overflow options is not read",
        ),
        (
            reading(
                &rolling(BlendRequest::new().edge(e.clone(), distances.clone())),
                Vary::Scale,
            ),
            "a reading of a round of two distances with rolling corners is not read",
        ),
        (
            reading(
                &rolling(BlendRequest::new().chain([e.clone()], law.clone())),
                Vary::Scale,
            ),
            "a reading of a round law with rolling corners is not read",
        ),
        (
            reading(
                &BlendRequest::new()
                    .chamfer([edges[0].clone()], 1.0)
                    .chamfer([edges[1].clone()], 2.0),
                Vary::Scale,
            ),
            "a scale reading of bevels of different measures is not read",
        ),
        (
            reading(&BlendRequest::new().round(edges.clone(), 0.0), Vary::Scale),
            "a scale reading of rounds of radius 0 reads nothing",
        ),
    ];
    for (got, text) in &cases {
        assert_eq!(got, &refused(text));
    }
}

/// `Blend` converts to the section it names, and the request's corners
/// default as the calls without them do.
#[test]
fn sections_and_defaults() {
    assert_eq!(ChamferCorner::default(), ChamferCorner::Facet);
    assert_eq!(Vary::default(), Vary::Scale);
    assert!(matches!(BlendSection::from(Blend::Round(2.0)), BlendSection::Round(r) if r == 2.0));
    assert!(matches!(
        BlendSection::from(Chamfer::Symmetric(1.0)),
        BlendSection::Chamfer(Chamfer::Symmetric(d)) if d == 1.0
    ));
}
