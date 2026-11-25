//! The largest measure several bevels take at once, confirmed against the
//! build ([`BlendSizing::Chamfers`]): at the reading the request builds
//! under its corner, and a [`LIMIT_TOLERANCE`] past it it is refused.
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{
    BlendLimit, BlendSizing, Chamfer, ChamferCorner, LIMIT_SAMPLES, LIMIT_TOLERANCE, Overflow,
    blend_limit, chamfer_edges_cornered, largest_chamfer_distance_all,
};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const CORNERS: [ChamferCorner; 2] = [ChamferCorner::Facet, ChamferCorner::Point];

fn boxed(dims: (f64, f64, f64)) -> (Model, Shape) {
    let mut m = Model::new();
    let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
    let s = ogeom_algo::make_box(&mut m, frame, dims, T).unwrap().shape;
    (m, s)
}

/// The edges of `s` every sample of which `on` takes.
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

fn symmetric(edges: &[Shape], d: f64) -> Vec<(Shape, Chamfer)> {
    edges
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(d)))
        .collect()
}

fn builds(m: &Model, s: &Shape, edges: &[Shape], d: f64, corner: ChamferCorner) -> bool {
    let mut scratch = m.clone();
    match chamfer_edges_cornered(&mut scratch, s, &symmetric(edges, d), corner, T) {
        Ok(built) => ogeom_algo::check(&scratch, &built.shape, T)
            .unwrap()
            .is_valid(),
        Err(_) => false,
    }
}

/// The confirmed reading of the symmetric bevels of `edges` under
/// `corner`, held to its promise: it builds valid and a
/// [`LIMIT_TOLERANCE`] past it is refused.
fn confirmed(m: &Model, s: &Shape, edges: &[Shape], corner: ChamferCorner) -> BlendLimit {
    let specs = symmetric(edges, 1.0);
    let found = blend_limit(
        m,
        s,
        &BlendSizing::Chamfers {
            specs: &specs,
            corner,
            overflow: Overflow::default(),
        },
        100.0,
        T,
    )
    .unwrap();
    let size = found.size.unwrap();
    assert!(builds(m, s, edges, size, corner), "{corner:?}: {found:?}");
    assert!(
        !builds(m, s, edges, size * (1.0 + LIMIT_TOLERANCE), corner),
        "{corner:?}: builds past the reading: {found:?}"
    );
    found
}

/// Both corners confirmed on `edges`, each read past `past` and within a
/// thousandth under `about`; [`largest_chamfer_distance_all`] answers the
/// facet's reading.
fn read_past_the_seat(m: &Model, s: &Shape, edges: &[Shape], past: f64, about: f64) {
    let mut facet = None;
    for corner in CORNERS {
        let found = confirmed(m, s, edges, corner);
        let size = found.size.unwrap();
        assert!(found.seat.unwrap() <= past, "{corner:?}: {found:?}");
        assert!(
            size > past * (1.0 + LIMIT_TOLERANCE) && (about * (1.0 - 1e-3)..=about).contains(&size),
            "{corner:?}: {found:?}"
        );
        if corner == ChamferCorner::Facet {
            facet = found.size;
        }
    }
    let all = largest_chamfer_distance_all(m, s, &symmetric(edges, 1.0), 100.0, T).unwrap();
    assert_eq!(all, facet);
}

/// The four top edges of a 20 x 20 x 10 box: the seat stops at the
/// sides' height, 10, but the bevels build on, notching over the sides,
/// until their legs on the top reach across it at 20.
#[test]
fn a_box_s_four_top_bevels_build_past_their_seat() {
    let (m, s) = boxed((20.0, 20.0, 10.0));
    let top = edges_where(&m, &s, |p| (p.z - 10.0).abs() < 1e-9);
    assert_eq!(top.len(), 4);
    read_past_the_seat(&m, &s, &top, 10.0, 20.0);
}

/// The top and bottom front edges of a 2 mm plate: past half the front
/// the two bevels meet along their ridge and build on, to about 21.
#[test]
fn two_bevels_across_a_plate_s_front_build_past_their_ridge() {
    let (m, s) = boxed((20.0, 20.0, 2.0));
    let front = edges_where(&m, &s, |p| {
        p.y.abs() < 1e-9 && (p.z.abs() < 1e-9 || (p.z - 2.0).abs() < 1e-9)
    });
    assert_eq!(front.len(), 2);
    read_past_the_seat(&m, &s, &front, 1.0, 21.0);
}

/// The two top edges of a 20 x 4 x 20 box, across its 4 mm top: past
/// half the top they meet along their ridge and build on, to about 22.
#[test]
fn two_bevels_across_a_narrow_top_build_past_their_ridge() {
    let (m, s) = boxed((20.0, 4.0, 20.0));
    let across = edges_where(&m, &s, |p| {
        (p.z - 20.0).abs() < 1e-9 && (p.y.abs() < 1e-9 || (p.y - 4.0).abs() < 1e-9)
    });
    assert_eq!(across.len(), 2);
    read_past_the_seat(&m, &s, &across, 2.0, 22.0);
}

/// Every edge of a 2 mm plate, its corners closed at a point: the seat,
/// half the thickness, is the limit, bracketed in two builds and the
/// samples.
#[test]
fn every_bevel_of_a_plate_meeting_at_points_stops_at_the_seat() {
    let (m, s) = boxed((20.0, 20.0, 2.0));
    let every = explore_unique(&m, &s, ShapeType::Edge).unwrap();
    let found = confirmed(&m, &s, &every, ChamferCorner::Point);
    let size = found.size.unwrap();
    assert!((1.0 - 1e-6..=1.0).contains(&size), "{found:?}");
    assert_eq!(found.probes, 2 + 2 * LIMIT_SAMPLES, "{found:?}");
}

#[test]
fn no_bevels_have_no_largest_measure() {
    let (m, s) = boxed((20.0, 20.0, 2.0));
    let sizing = BlendSizing::Chamfers {
        specs: &[],
        corner: ChamferCorner::Facet,
        overflow: Overflow::default(),
    };
    assert!(blend_limit(&m, &s, &sizing, 10.0, T).is_err());
}
