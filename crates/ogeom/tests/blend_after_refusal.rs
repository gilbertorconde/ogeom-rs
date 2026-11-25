//! A blend call on a model that holds other shapes, or what earlier calls
//! left, beside the given solid answers as on a model holding only that
//! solid, read back from its native text: the same volume, or the same
//! refusal, never a panic.
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::Built;
use ogeom::core::{OgeomResult, Tolerances};
use ogeom::fillet::{
    Chamfer, ChamferCorner, FilletCorner, Overflow, chamfer_edges_overflow, fillet_edges,
    fillet_edges_corner, fillet_edges_overflow,
};
use ogeom::geom::Curve3d as _;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// A 40 by 40 by 4 plate with a boss of radius 5, 5 tall, and a 0.5 square
/// pocket cut down through the boss's top just inside its rim, built in a
/// model that keeps the plate, the boss and the pocket.
fn bossed() -> (Model, Shape) {
    let mut m = Model::new();
    let plate = ogeom::algo::make_box(&mut m, at(-20.0, -20.0, 0.0), (40.0, 40.0, 4.0), T)
        .unwrap()
        .shape;
    let boss = ogeom::algo::make_cylinder(&mut m, at(0.0, 0.0, 4.0), 5.0, 5.0, T)
        .unwrap()
        .shape;
    let s = ogeom::boolean::fuse(&mut m, &plate, &boss, T)
        .unwrap()
        .shape;
    let pocket = ogeom::algo::make_box(&mut m, at(4.3, -0.25, 6.0), (0.5, 0.5, 4.0), T)
        .unwrap()
        .shape;
    let s = ogeom::boolean::cut(&mut m, &s, &pocket, T).unwrap().shape;
    (m, s)
}

/// The model holding `s` alone, read back from its native text.
fn alone(m: &Model, s: &Shape) -> (Model, Shape) {
    let text = ogeom::io::native::write(
        m,
        std::slice::from_ref(s),
        ogeom::io::native::WriteOptions::default(),
    )
    .unwrap();
    let (m, roots) = ogeom::io::native::read(&text).unwrap();
    (m, roots[0].clone())
}

/// The edges of `s` every sample of which `on` takes, in the order of
/// their middles.
fn edges_where(m: &Model, s: &Shape, on: impl Fn(Point) -> bool) -> Vec<Shape> {
    let mut found: Vec<(Shape, [f64; 3])> = explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter_map(|e| {
            let (c, (a, b)) = ogeom::algo::placed_edge_curve(m, &e, T).unwrap().unwrap();
            let ps: Vec<Point> = [0.0, 0.3, 0.5, 0.7, 1.0]
                .iter()
                .map(|k| c.point_at(a + (b - a) * k, T).unwrap())
                .collect();
            ps.iter()
                .all(|p| on(*p))
                .then(|| (e, [ps[2].x, ps[2].y, ps[2].z]))
        })
        .collect();
    found.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    found.into_iter().map(|(e, _)| e).collect()
}

/// The two half circles of the boss's top rim.
fn rim(m: &Model, s: &Shape) -> Vec<Shape> {
    let rim = edges_where(m, s, |p| {
        (p.z - 9.0).abs() < 1e-6 && (p.x.hypot(p.y) - 5.0).abs() < 1e-6
    });
    assert_eq!(rim.len(), 2);
    rim
}

/// The plate's top edge along `y = -20`.
fn plate_edge(m: &Model, s: &Shape) -> Vec<Shape> {
    let e = edges_where(m, s, |p| {
        (p.z - 4.0).abs() < 1e-6 && (p.y + 20.0).abs() < 1e-6
    });
    assert_eq!(e.len(), 1);
    e
}

fn no_trim() -> Overflow {
    let mut o = Overflow::default();
    o.trim_through = false;
    o
}

/// What a call answered: the result's volume, or the refusal's text.
#[derive(Debug)]
enum Answer {
    Built(f64),
    Refused(String),
}

fn answer(m: &Model, r: OgeomResult<Built>) -> Answer {
    match r {
        Ok(b) => Answer::Built(
            ogeom::algo::volume_properties(m, &b.shape, Deflection::default(), T)
                .unwrap()
                .mass,
        ),
        Err(e) => Answer::Refused(e.to_string()),
    }
}

fn assert_same(got: &Answer, want: &Answer, what: &str) {
    let same = match (got, want) {
        (Answer::Built(x), Answer::Built(y)) => (x - y).abs() <= 1e-6 * y.abs(),
        (Answer::Refused(x), Answer::Refused(y)) => x == y,
        _ => false,
    };
    assert!(same, "{what}: {got:?} here, {want:?} on the solid alone");
}

type Call = fn(&mut Model, &Shape, &[Shape]) -> OgeomResult<Built>;

fn overflow_default(m: &mut Model, s: &Shape, e: &[Shape]) -> OgeomResult<Built> {
    fillet_edges_overflow(m, s, e, 1.0, Overflow::default(), T)
}

fn overflow_no_trim(m: &mut Model, s: &Shape, e: &[Shape]) -> OgeomResult<Built> {
    fillet_edges_overflow(m, s, e, 1.0, no_trim(), T)
}

fn plain(m: &mut Model, s: &Shape, e: &[Shape]) -> OgeomResult<Built> {
    fillet_edges(m, s, e, 1.0, T)
}

fn setback(m: &mut Model, s: &Shape, e: &[Shape]) -> OgeomResult<Built> {
    fillet_edges_corner(m, s, e, 1.0, FilletCorner::Setback, T)
}

fn rolling(m: &mut Model, s: &Shape, e: &[Shape]) -> OgeomResult<Built> {
    fillet_edges_corner(m, s, e, 1.0, FilletCorner::Rolling, T)
}

fn chamfer(m: &mut Model, s: &Shape, e: &[Shape]) -> OgeomResult<Built> {
    let specs: Vec<(Shape, Chamfer)> = e
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(0.5)))
        .collect();
    chamfer_edges_overflow(m, s, &specs, ChamferCorner::Point, Overflow::default(), T)
}

/// The rim refused whole with trim through off, then one half of it on a
/// copy of that model, answers as on the solid alone.
fn one_half_after_the_rim_is_refused(call: Call) {
    let (mut m, s) = bossed();
    let (mut lone, lone_s) = alone(&m, &s);
    let lone_rim = rim(&lone, &lone_s);
    let want = call(&mut lone, &lone_s, &lone_rim[..1]);
    let want = answer(&lone, want);

    let halves = rim(&m, &s);
    assert!(overflow_no_trim(&mut m, &s, &halves).is_err());
    let mut again = m.clone();
    let got = call(&mut again, &s, &halves[..1]);
    assert_same(&answer(&again, got), &want, "one half");
}

#[test]
fn overflow_after_refusal() {
    one_half_after_the_rim_is_refused(overflow_default);
}

#[test]
fn overflow_without_trim_after_refusal() {
    one_half_after_the_rim_is_refused(overflow_no_trim);
}

#[test]
fn fillet_after_refusal() {
    one_half_after_the_rim_is_refused(plain);
}

#[test]
fn setback_corner_after_refusal() {
    one_half_after_the_rim_is_refused(setback);
}

#[test]
fn rolling_corner_after_refusal() {
    one_half_after_the_rim_is_refused(rolling);
}

#[test]
fn chamfer_after_refusal() {
    one_half_after_the_rim_is_refused(chamfer);
}

/// Refused requests and valid ones in turn on one model, each valid one
/// on the original solid answering as on that solid alone.
#[test]
#[ignore = "heavy: eight blend calls, each against the solid alone"]
fn valid_calls_after_refusals_on_one_model() {
    let (mut m, s) = bossed();
    let (lone, lone_s) = alone(&m, &s);
    let pick = |m: &Model, s: &Shape, which: usize| -> Vec<Shape> {
        match which {
            0 => rim(m, s)[..1].to_vec(),
            1 => rim(m, s)[1..].to_vec(),
            _ => plate_edge(m, s),
        }
    };
    let steps: [(&str, Call, usize); 7] = [
        ("one half, overflow", overflow_default, 0),
        ("other half, fillet", plain, 1),
        ("one half, chamfer", chamfer, 0),
        ("plate edge, fillet", plain, 2),
        ("other half, rolling", rolling, 1),
        ("plate edge, chamfer", chamfer, 2),
        ("one half, setback", setback, 0),
    ];
    let whole = rim(&m, &s);
    assert!(overflow_no_trim(&mut m, &s, &whole).is_err());
    for (k, (what, call, which)) in steps.iter().enumerate() {
        if k == 3 {
            assert!(overflow_no_trim(&mut m, &s, &whole).is_err());
        }
        let mut fresh = lone.clone();
        let asked = pick(&fresh, &lone_s, *which);
        let want = call(&mut fresh, &lone_s, &asked);
        let want = answer(&fresh, want);
        let edges = pick(&m, &s, *which);
        let got = call(&mut m, &s, &edges);
        assert_same(&answer(&m, got), &want, what);
    }
}

/// The boss's top face, the plane beside the rim.
fn top(m: &Model, s: &Shape, edge: &Shape) -> Shape {
    ogeom::topo::ancestors_of(m, s, edge, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let d = m.node(f).unwrap().data().as_face().unwrap();
            matches!(
                m.geometry().surface(d.surface),
                Some(ogeom::geom::SurfaceGeometry::Plane(_))
            )
        })
        .unwrap()
}

/// Each half of the rim rounded alone by `fillet_edge`, and by
/// `fillet_edge_distances` at equal distances, with no request before it:
/// as on the solid alone.
#[test]
fn each_half_rounded_alone_as_on_the_solid_alone() {
    let (m, s) = bossed();
    let (lone, lone_s) = alone(&m, &s);
    type One = fn(&mut Model, &Shape, &Shape) -> OgeomResult<Built>;
    let calls: [(&str, One); 2] = [
        ("fillet_edge", |m, s, e| {
            ogeom::fillet::fillet_edge(m, s, e, 1.0, T)
        }),
        ("equal distances", |m, s, e| {
            let face = top(m, s, e);
            ogeom::fillet::fillet_edge_distances(
                m,
                s,
                e,
                &face,
                1.0,
                1.0,
                ogeom::fillet::RoundProfile::Elliptical,
                T,
            )
        }),
    ];
    for (what, call) in calls {
        for k in 0..2 {
            let mut fresh = lone.clone();
            let edge = rim(&fresh, &lone_s)[k].clone();
            let want = call(&mut fresh, &lone_s, &edge);
            let want = answer(&fresh, want);
            let mut here = m.clone();
            let edge = rim(&here, &s)[k].clone();
            let got = call(&mut here, &s, &edge);
            assert_same(&answer(&here, got), &want, &format!("{what}, half {k}"));
        }
    }
}
