//! Variable-radius rounds under a radius law: the straight planar seat's
//! exact envelope, and closed curved seats marched with the law's radius.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{
    RadiusLaw, analyse_blend, fillet_edge, fillet_edge_variable, fillet_edges_variable,
};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Orientation, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// The edge of `s` passing nearest `p`, 64 samples along each.
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

fn volume(m: &Model, s: &Shape, chord: f64) -> f64 {
    let deflection = ogeom_mesh::Deflection {
        chord,
        ..ogeom_mesh::Deflection::default()
    };
    ogeom_algo::volume_properties(m, s, deflection, T)
        .unwrap()
        .mass
}

fn valid(m: &Model, s: &Shape) {
    let diagnosis = ogeom_algo::check(m, s, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
}

/// The law `(0, 1), (0.5, 2), (1, 1)` in closed form: on an open chain
/// the end slopes are the secants (2 and -2), the middle one level where
/// the secants change sign, so each half is the Hermite cubic between.
fn hump_open(s: f64) -> f64 {
    let (t, m0, m1, y0, y1) = if s <= 0.5 {
        (s / 0.5, 2.0, 0.0, 1.0, 2.0)
    } else {
        ((s - 0.5) / 0.5, 0.0, -2.0, 2.0, 1.0)
    };
    let h = 0.5;
    let (t2, t3) = (t * t, t * t * t);
    (2.0 * t3 - 3.0 * t2 + 1.0) * y0
        + (t3 - 2.0 * t2 + t) * m0 * h
        + (3.0 * t2 - 2.0 * t3) * y1
        + (t3 - t2) * m1 * h
}

/// The same points on a closed chain: 0 and 1 are one point, and every
/// slope is level, the secants either side of each point changing sign.
fn hump_closed(s: f64) -> f64 {
    let t = if s <= 0.5 { s / 0.5 } else { (s - 0.5) / 0.5 };
    let (y0, y1) = if s <= 0.5 { (1.0, 2.0) } else { (2.0, 1.0) };
    let (t2, t3) = (t * t, t * t * t);
    (2.0 * t3 - 3.0 * t2 + 1.0) * y0 + (3.0 * t2 - 2.0 * t3) * y1
}

fn hump() -> RadiusLaw {
    RadiusLaw::Points(vec![(0.0, 1.0), (0.5, 2.0), (1.0, 1.0)])
}

#[test]
fn a_linear_law_on_a_straight_edge_is_the_linear_fillet() {
    let mut model = Model::new();
    let block = ogeom_algo::make_box(&mut model, at(0.0, 0.0, 0.0), (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edge = edge_near(&model, &block, [10.0, 0.0, 10.0]);
    let law = RadiusLaw::Linear {
        start: 1.0,
        end: 2.0,
    };
    let by_law = fillet_edges_variable(&mut model, &block, std::slice::from_ref(&edge), &law, T)
        .unwrap()
        .shape;
    // The chain starts where the edge, as oriented, does.
    let (from, to) = if edge.orientation() == Orientation::Reversed {
        (2.0, 1.0)
    } else {
        (1.0, 2.0)
    };
    let linear = fillet_edge_variable(&mut model, &block, &edge, from, to, T)
        .unwrap()
        .shape;
    valid(&model, &by_law);
    let (a, b) = (volume(&model, &by_law, 1e-3), volume(&model, &linear, 1e-3));
    assert!((a - b).abs() < 1e-12, "{a} against {b}");
}

#[test]
fn a_law_of_points_on_a_straight_edge_is_its_exact_envelope() {
    // A 20 x 20 x 10 block, its top front edge from x = 0 to x = 20,
    // rounded by the hump law: the removed sliver is (1 - pi/4) times the
    // integral of the radius squared along the edge, and each section's
    // arc touches both faces at the law's radius from the edge.
    let mut model = Model::new();
    let block = ogeom_algo::make_box(&mut model, at(0.0, 0.0, 0.0), (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edge = edge_near(&model, &block, [10.0, 0.0, 10.0]);
    let rounded =
        fillet_edges_variable(&mut model, &block, std::slice::from_ref(&edge), &hump(), T).unwrap();
    valid(&model, &rounded.shape);
    assert!(rounded.history.is_deleted(&edge));

    // The chain's fraction at x, by the edge's orientation.
    let reversed = {
        let (c, (lo, _)) = ogeom_algo::placed_edge_curve(&model, &edge, T)
            .unwrap()
            .unwrap();
        let start_x = c.point_at(lo, T).unwrap().x;
        (start_x > 10.0) != (edge.orientation() == Orientation::Reversed)
    };
    let fraction = |x: f64| if reversed { 1.0 - x / 20.0 } else { x / 20.0 };

    let n = 2000;
    let mut integral = 0.0;
    for i in 0..n {
        let s = (f64::from(i) + 0.5) / f64::from(n);
        integral += hump_open(s).powi(2) / f64::from(n);
    }
    let exact = 4000.0 - 20.0 * (1.0 - core::f64::consts::PI / 4.0) * integral;
    let measured = volume(&model, &rounded.shape, 1e-4);
    assert!(
        (measured - exact).abs() < 2e-3,
        "{measured} against {exact}"
    );

    // The blend face: the one B-spline face, its section at each x the
    // quarter arc about (x, r, 10 - r) of radius r.
    let blend = explore_unique(&model, &rounded.shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            model
                .node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::BSpline(_)))
        })
        .expect("the blend face");
    let surface = model
        .node(&blend)
        .and_then(|n| n.data().as_face())
        .and_then(|d| model.geometry().surface(d.surface))
        .cloned()
        .unwrap();
    for i in 0..=16 {
        let x = 20.0 * f64::from(i) / 16.0;
        let r = hump_open(fraction(x));
        for k in 0..=4 {
            let a = core::f64::consts::FRAC_PI_2 * f64::from(k) / 4.0;
            let p = Point::new(x, r - r * a.cos(), 10.0 - r + r * a.sin());
            let off = ogeom_algo::project_on_surface(&surface, p, 24, T)
                .unwrap()
                .point
                .distance(p);
            assert!(off < 1e-9, "{p:?} is {off} off the blend");
        }
    }
    for contact in analyse_blend(&model, &rounded.shape, &blend, 33, T).unwrap() {
        if contact.tangency_error < 1.0 {
            assert!(contact.tangency_error < 1e-9, "{contact:?}");
        }
    }
}

/// A contact's fraction along the chain and its distance from the rim.
type Setback = (f64, f64);

/// The two rails of a rim's round on a drum of radius `drum` and height
/// `height`: the points on edges of `shape` on the top inside the rim,
/// and on the wall below it, each with the chain's fraction it stands at.
fn rim_contacts(
    model: &Model,
    shape: &Shape,
    drum: f64,
    height: f64,
    fraction: impl Fn(Point) -> f64,
) -> (Vec<Setback>, Vec<Setback>) {
    let (mut top, mut wall) = (Vec::new(), Vec::new());
    for e in explore_unique(model, shape, ShapeType::Edge).unwrap() {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(model, &e, T).unwrap() else {
            continue;
        };
        // A rail lies wholly on the top inside the rim, or wholly on the
        // wall below it.
        let points: Vec<Point> = (0..=64)
            .map(|i| c.point_at(a + (b - a) * f64::from(i) / 64.0, T).unwrap())
            .collect();
        let on_top = |p: &Point| {
            let rho = p.x.hypot(p.y);
            (p.z - height).abs() < 1e-3 && rho < drum - 1e-3 && rho > drum - 3.0
        };
        let on_wall = |p: &Point| {
            (p.x.hypot(p.y) - drum).abs() < 1e-3 && p.z < height - 1e-3 && p.z > height - 3.0
        };
        if points.iter().all(on_top) {
            top.extend(points.iter().map(|p| (fraction(*p), drum - p.x.hypot(p.y))));
        } else if points.iter().all(on_wall) {
            wall.extend(points.iter().map(|p| (fraction(*p), height - p.z)));
        }
    }
    (top, wall)
}

#[test]
fn a_law_of_points_rounds_a_drum_s_rim() {
    let mut model = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, at(0.0, 0.0, 0.0), 10.0, 20.0, T)
        .unwrap()
        .shape;
    let rim = edge_near(&model, &drum, [0.0, -10.0, 20.0]);
    let rounded =
        fillet_edges_variable(&mut model, &drum, std::slice::from_ref(&rim), &hump(), T).unwrap();
    valid(&model, &rounded.shape);
    assert!(rounded.history.is_deleted(&rim));

    // The chain's fraction at a point: its angle round the axis from the
    // rim's start, the way the rim runs.
    let (start, ahead) = {
        let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(&model, &rim, T)
            .unwrap()
            .unwrap();
        let (t, sense) = if rim.orientation() == Orientation::Reversed {
            (hi, -1.0)
        } else {
            (lo, 1.0)
        };
        (c.point_at(t, T).unwrap(), c.d1_at(t, T).unwrap() * sense)
    };
    let fraction = |p: Point| {
        let a0 = start.y.atan2(start.x);
        let turn = if (start.x * ahead.y - start.y * ahead.x) > 0.0 {
            1.0
        } else {
            -1.0
        };
        let a = p.y.atan2(p.x);
        ((a - a0) * turn).rem_euclid(core::f64::consts::TAU) / core::f64::consts::TAU
    };
    let (top, wall) = rim_contacts(&model, &rounded.shape, 10.0, 20.0, fraction);
    assert!(
        top.len() > 16 && wall.len() > 16,
        "{} and {}",
        top.len(),
        wall.len()
    );
    for (s, setback) in top.iter().chain(&wall) {
        let law = hump_closed(*s);
        assert!(
            (setback - law).abs() < 1e-3,
            "at {s} the contact stands {setback} from the rim, the law {law}"
        );
    }

    // Between the constant rims of the law's least and largest radii.
    let mut constant = |r: f64| {
        let d = ogeom_algo::make_cylinder(&mut model, at(0.0, 0.0, 0.0), 10.0, 20.0, T)
            .unwrap()
            .shape;
        let e = edge_near(&model, &d, [0.0, -10.0, 20.0]);
        fillet_edge(&mut model, &d, &e, r, T).unwrap().shape
    };
    let (one, two) = (constant(1.0), constant(2.0));
    let v = volume(&model, &rounded.shape, 1e-3);
    let (v1, v2) = (volume(&model, &one, 1e-3), volume(&model, &two, 1e-3));
    assert!(v2 < v && v < v1, "{v2} < {v} < {v1}");
}

#[test]
fn a_law_of_points_rounds_a_marched_seat_tangent_to_both_faces() {
    // A drum cut by a plane leaning twenty degrees: its top's rim is an
    // oblique ellipse, closed, which only the march rounds.
    let mut model = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, at(0.0, 0.0, 0.0), 5.0, 20.0, T)
        .unwrap()
        .shape;
    let lean = 20.0_f64.to_radians();
    let normal = ogeom_math::Vector::new(lean.sin(), 0.0, lean.cos());
    let across = ogeom_math::Vector::new(lean.cos(), 0.0, -lean.sin());
    let side = normal.cross(across);
    let corner = Point::new(0.0, 0.0, 15.0) - across * 20.0 - side * 20.0;
    let frame = Frame::new(
        corner,
        Direction::new(normal, T).unwrap(),
        Direction::new(across, T).unwrap(),
        T,
    )
    .unwrap();
    let above = ogeom_algo::make_box(&mut model, frame, (40.0, 40.0, 20.0), T)
        .unwrap()
        .shape;
    let joined = ogeom_bool::cut(&mut model, &drum, &above, T).unwrap().shape;
    // The rim: the two arcs the boolean split it into, end to end.
    let arcs: Vec<Shape> = explore_unique(&model, &joined, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(&model, e, T)
                .unwrap()
                .unwrap();
            c.point_at(lo, T).unwrap().z > 10.0 && c.point_at(hi, T).unwrap().z > 10.0
        })
        .collect();
    assert_eq!(arcs.len(), 2, "the oblique rim is two arcs");
    let ends = |e: &Shape| ogeom_algo::edge_vertices(&model, e).unwrap().unwrap();
    let chain = if ends(&arcs[0]).1.is_same(&ends(&arcs[1]).0) {
        vec![arcs[0].clone(), arcs[1].clone()]
    } else {
        vec![arcs[0].clone(), arcs[1].reversed()]
    };
    let seat = chain[0].clone();

    let law = RadiusLaw::Points(vec![(0.0, 0.5), (0.5, 1.0)]);
    // One arc alone is a crease that closes through the other: its law
    // would be read round the whole loop the march walks.
    match fillet_edges_variable(&mut model, &joined, &chain[..1], &law, T) {
        Err(e) => assert!(e.to_string().contains("closes through"), "{e}"),
        Ok(_) => panic!("one arc of a closed crease took a law of its own"),
    }
    let rounded = fillet_edges_variable(&mut model, &joined, &chain, &law, T).unwrap();
    valid(&model, &rounded.shape);
    assert!(rounded.history.is_deleted(&seat));

    let blends: Vec<Shape> = rounded
        .history
        .generated(&seat)
        .iter()
        .filter(|f| model.kind_of(f).is_ok_and(|k| k == ShapeType::Face))
        .cloned()
        .collect();
    assert_eq!(blends.len(), 1, "one band");
    // Tangent to the top and the wall along every piece of its rails, as
    // closely as the constant round's fitted band is.
    let worst = |model: &Model, shape: &Shape, band: &Shape| -> (usize, f64) {
        let contacts = analyse_blend(model, shape, band, 33, T).unwrap();
        let meets: Vec<_> = contacts.iter().filter(|c| c.tangency_error < 1.0).collect();
        let mut faces: Vec<Shape> = Vec::new();
        for c in &meets {
            if !faces.iter().any(|f| f.is_same(&c.neighbour)) {
                faces.push(c.neighbour.clone());
            }
        }
        let error = meets.iter().map(|c| c.tangency_error).fold(0.0, f64::max);
        (faces.len(), error)
    };
    let (faces, error) = worst(&model, &rounded.shape, &blends[0]);
    assert_eq!(faces, 2, "the band meets the top and the wall");

    // Taken off, and less than the round at the law's largest radius
    // takes.
    let before = volume(&model, &joined, 1e-3);
    let removed = before - volume(&model, &rounded.shape, 1e-3);
    let largest = fillet_edge(&mut model, &joined, &seat, 1.0, T).unwrap();
    let band: Vec<Shape> = largest
        .history
        .generated(&seat)
        .iter()
        .filter(|f| model.kind_of(f).is_ok_and(|k| k == ShapeType::Face))
        .cloned()
        .collect();
    let (_, constant) = worst(&model, &largest.shape, &band[0]);
    assert!(
        error <= constant * 1.001,
        "tangent within {error} rad, the constant band within {constant}"
    );
    let most = before - volume(&model, &largest.shape, 1e-3);
    assert!(
        removed > 0.0 && removed < most,
        "removed {removed} of at most {most}"
    );
}

/// A drum of radius 5 cut in half through its axis: its top's rim is an
/// open half circle, a crease that ends where the cut meets it, so the
/// round runs out through the cut at both ends.
fn half_drum(model: &mut Model) -> (Shape, Shape) {
    let drum = ogeom_algo::make_cylinder(model, at(0.0, 0.0, 0.0), 5.0, 10.0, T)
        .unwrap()
        .shape;
    let away = ogeom_algo::make_box(model, at(-10.0, -10.0, -1.0), (20.0, 10.0, 12.0), T)
        .unwrap()
        .shape;
    let half = ogeom_bool::cut(model, &drum, &away, T).unwrap().shape;
    let rim = edge_near(model, &half, [0.0, 5.0, 10.0]);
    (half, rim)
}

/// An open curved edge under a law of points: the law is laid along the
/// edge's own arc length and held at its end radius past each end, where
/// the ball runs out at that one radius, as a constant round's does. Every
/// contact on both rails stands at the law's setback for its place along
/// the arc, and the round removes more than the constant round of the
/// law's least radius and less than its largest's.
#[test]
fn a_law_of_points_rounds_an_open_arc_running_out_at_both_ends() {
    let mut model = Model::new();
    let (half, rim) = half_drum(&mut model);
    let rounded =
        fillet_edges_variable(&mut model, &half, std::slice::from_ref(&rim), &hump(), T).unwrap();
    valid(&model, &rounded.shape);
    assert!(rounded.history.is_deleted(&rim));

    // The arc's fraction at a point: its angle round the axis from the
    // rim's start, the way the rim runs, over the half turn.
    let (start, ahead) = {
        let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(&model, &rim, T)
            .unwrap()
            .unwrap();
        let (t, sense) = if rim.orientation() == Orientation::Reversed {
            (hi, -1.0)
        } else {
            (lo, 1.0)
        };
        (c.point_at(t, T).unwrap(), c.d1_at(t, T).unwrap() * sense)
    };
    let fraction = |p: Point| {
        let a0 = start.y.atan2(start.x);
        let turn = if (start.x * ahead.y - start.y * ahead.x) > 0.0 {
            1.0
        } else {
            -1.0
        };
        let a = p.y.atan2(p.x);
        ((a - a0) * turn).rem_euclid(core::f64::consts::TAU) / core::f64::consts::PI
    };
    let (top, wall) = rim_contacts(&model, &rounded.shape, 5.0, 10.0, fraction);
    assert!(
        top.len() > 16 && wall.len() > 16,
        "{} and {}",
        top.len(),
        wall.len()
    );
    for (s, setback) in top.iter().chain(&wall) {
        let law = hump_open(s.clamp(0.0, 1.0));
        assert!(
            (setback - law).abs() < 1e-3,
            "at {s} the contact stands {setback} from the rim, the law {law}"
        );
    }

    let mut constant = |r: f64| {
        let (d, e) = half_drum(&mut model);
        fillet_edge(&mut model, &d, &e, r, T).unwrap().shape
    };
    let (one, two) = (constant(1.0), constant(2.0));
    let v = volume(&model, &rounded.shape, 1e-3);
    let (v1, v2) = (volume(&model, &one, 1e-3), volume(&model, &two, 1e-3));
    assert!(v2 < v && v < v1, "{v2} < {v} < {v1}");
}

/// An open chain of arcs of one curve, a crease split at a vertex, under a
/// law of points: the half of a radius-5 drum across its seam, whose top
/// rim is two arcs of one circle meeting at the seam's vertex. It is
/// rounded as the one seat it is, the law laid along both arcs' length
/// and the ball running out through the cut at both ends.
#[test]
fn a_law_of_points_rounds_an_open_chain_of_arcs_of_one_circle() {
    let mut model = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, at(0.0, 0.0, 0.0), 5.0, 10.0, T)
        .unwrap()
        .shape;
    let away = ogeom_algo::make_box(&mut model, at(-10.0, -10.0, -1.0), (10.0, 20.0, 12.0), T)
        .unwrap()
        .shape;
    let half = ogeom_bool::cut(&mut model, &drum, &away, T).unwrap().shape;
    let arcs: Vec<Shape> = explore_unique(&model, &half, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(&model, e, T)
                .unwrap()
                .unwrap();
            let mid = c.point_at(f64::midpoint(lo, hi), T).unwrap();
            (mid.z - 10.0).abs() < 1e-6 && (mid.x.hypot(mid.y) - 5.0).abs() < 1e-6
        })
        .collect();
    assert_eq!(arcs.len(), 2, "the rim is two arcs, split at the seam");
    let ends = |e: &Shape| ogeom_algo::edge_vertices(&model, e).unwrap().unwrap();
    let chain = if ends(&arcs[0]).1.is_same(&ends(&arcs[1]).0) {
        vec![arcs[0].clone(), arcs[1].clone()]
    } else {
        vec![arcs[1].clone(), arcs[0].clone()]
    };
    let rounded = fillet_edges_variable(&mut model, &half, &chain, &hump(), T).unwrap();
    valid(&model, &rounded.shape);
    assert!(chain.iter().all(|e| rounded.history.is_deleted(e)));

    // The chain's fraction at a point: its angle round the axis from the
    // chain's start, the way it runs, over the half turn.
    let (start, ahead) = {
        let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(&model, &chain[0], T)
            .unwrap()
            .unwrap();
        let (t, sense) = if chain[0].orientation() == Orientation::Reversed {
            (hi, -1.0)
        } else {
            (lo, 1.0)
        };
        (c.point_at(t, T).unwrap(), c.d1_at(t, T).unwrap() * sense)
    };
    let fraction = |p: Point| {
        let a0 = start.y.atan2(start.x);
        let turn = if (start.x * ahead.y - start.y * ahead.x) > 0.0 {
            1.0
        } else {
            -1.0
        };
        let a = p.y.atan2(p.x);
        ((a - a0) * turn).rem_euclid(core::f64::consts::TAU) / core::f64::consts::PI
    };
    let (top, wall) = rim_contacts(&model, &rounded.shape, 5.0, 10.0, fraction);
    assert!(
        top.len() > 16 && wall.len() > 16,
        "{} and {}",
        top.len(),
        wall.len()
    );
    for (s, setback) in top.iter().chain(&wall) {
        let law = hump_open(s.clamp(0.0, 1.0));
        assert!(
            (setback - law).abs() < 1e-3,
            "at {s} the contact stands {setback} from the rim, the law {law}"
        );
    }
}

/// A straight edge along a drum's ruling under a law of points: a drum of
/// radius 5, 10 tall, cut flat at x = 3, its flat meeting it along a ruling
/// at y = 4. The ball of radius `r` seated on the flat and inside the drum
/// has its centre at x = 3 - r, 5 - r from the axis, so its contact on the
/// flat stands at y = sqrt((5 - r)^2 - (3 - r)^2) and on the drum at
/// x = 5 (3 - r) / (5 - r); every contact on both rails is held there at
/// the law's radius for its height.
#[test]
fn a_law_of_points_rounds_a_flat_meeting_a_drum_along_its_ruling() {
    let mut model = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, at(0.0, 0.0, 0.0), 5.0, 10.0, T)
        .unwrap()
        .shape;
    let away = ogeom_algo::make_box(&mut model, at(3.0, -10.0, -1.0), (10.0, 20.0, 12.0), T)
        .unwrap()
        .shape;
    let flat = ogeom_bool::cut(&mut model, &drum, &away, T).unwrap().shape;
    let ruling = edge_near(&model, &flat, [3.0, 4.0, 5.0]);
    let rounded =
        fillet_edges_variable(&mut model, &flat, std::slice::from_ref(&ruling), &hump(), T)
            .unwrap();
    valid(&model, &rounded.shape);
    assert!(rounded.history.is_deleted(&ruling));

    // The fraction along the ruling at a height, the way it runs.
    let upward = {
        let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(&model, &ruling, T)
            .unwrap()
            .unwrap();
        let rises = c.point_at(hi, T).unwrap().z > c.point_at(lo, T).unwrap().z;
        rises != (ruling.orientation() == Orientation::Reversed)
    };
    let fraction = |z: f64| if upward { z / 10.0 } else { 1.0 - z / 10.0 };
    let (mut on_flat, mut on_drum) = (0, 0);
    for e in explore_unique(&model, &rounded.shape, ShapeType::Edge).unwrap() {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(&model, &e, T).unwrap() else {
            continue;
        };
        let points: Vec<Point> = (0..=32)
            .map(|i| c.point_at(a + (b - a) * f64::from(i) / 32.0, T).unwrap())
            .collect();
        // Rails run up the ruling's side from end face to end face: on the
        // flat beside it, or on the drum beside it (a fitted curve, held to
        // the drum within the band's fit).
        let (low, high) = points
            .iter()
            .fold((f64::MAX, f64::MIN), |(l, h), p| (l.min(p.z), h.max(p.z)));
        if high - low < 5.0 || points.iter().any(|p| p.y < 0.5) {
            continue;
        }
        if points
            .iter()
            .all(|p| (p.x - 3.0).abs() < 1e-6 && p.y < 4.0 - 1e-3)
        {
            for p in &points {
                let r = hump_open(fraction(p.z));
                let want = ((5.0 - r).powi(2) - (3.0 - r).powi(2)).sqrt();
                assert!(
                    (p.y - want).abs() < 1e-3,
                    "on the flat at z {}: y {} against {want}",
                    p.z,
                    p.y
                );
                on_flat += 1;
            }
        } else if points
            .iter()
            .all(|p| (p.x.hypot(p.y) - 5.0).abs() < 1e-3 && p.x < 3.0 - 1e-3)
        {
            for p in &points {
                let r = hump_open(fraction(p.z));
                let want = 5.0 * (3.0 - r) / (5.0 - r);
                assert!(
                    (p.x - want).abs() < 1e-3,
                    "on the drum at z {}: x {} against {want}",
                    p.z,
                    p.x
                );
                on_drum += 1;
            }
        }
    }
    assert!(on_flat > 16 && on_drum > 16, "{on_flat} and {on_drum}");
}

/// Two drums of radius 5 crossing square, the upright one about z and the
/// other about the line y = 0, z = 10, fused, and an edge of their seam:
/// an arc of an ellipse `z = 10 ± x` on the upright drum, ending at a pole
/// (0, ±5, 10) where the drums turn tangent.
///
/// The ball rolls outside both drums, its centre 5 + r from both axes; a
/// touch on the upright drum at (x, y, z) then puts it at radius
/// `5 |z - 10| / |x| - 5`, a touch on the other at `5 |x| / |z - 10| - 5`.
/// Each station's section stands square to the seam at its first touch,
/// so every rail point is held to the law at the seam's arc fraction of
/// its own foot there or of its ball's other touch's.
#[test]
fn a_law_of_points_rounds_a_seam_pinched_where_two_drums_turn_tangent() {
    let mut model = Model::new();
    let upright = ogeom_algo::make_cylinder(&mut model, at(0.0, 0.0, 0.0), 5.0, 20.0, T)
        .unwrap()
        .shape;
    let across_frame =
        Frame::new(Point::new(-20.0, 0.0, 10.0), Direction::X, Direction::Y, T).unwrap();
    let across = ogeom_algo::make_cylinder(&mut model, across_frame, 5.0, 40.0, T)
        .unwrap()
        .shape;
    let joined = ogeom_bool::fuse(&mut model, &upright, &across, T)
        .unwrap()
        .shape;
    let seam = explore_unique(&model, &joined, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(&model, e, T)
                .unwrap()
                .unwrap();
            [a, b].iter().any(|t| {
                let p = c.point_at(*t, T).unwrap();
                p.x.abs() < 1e-6 && (p.z - 10.0).abs() < 1e-6
            }) && (0..=8).all(|i| {
                let p = c.point_at(a + (b - a) * f64::from(i) / 8.0, T).unwrap();
                (p.x.hypot(p.y) - 5.0).abs() < 1e-6 && (p.y.hypot(p.z - 10.0) - 5.0).abs() < 1e-6
            })
        })
        .unwrap();

    // The seam sampled by arc length, the way the edge runs.
    let samples: Vec<(Point, f64)> = {
        let (c, (a, b)) = ogeom_algo::placed_edge_curve(&model, &seam, T)
            .unwrap()
            .unwrap();
        let mut points: Vec<Point> = (0..=20_000)
            .map(|i| {
                c.point_at(a + (b - a) * f64::from(i) / 20_000.0, T)
                    .unwrap()
            })
            .collect();
        if seam.orientation() == Orientation::Reversed {
            points.reverse();
        }
        let mut run = vec![0.0];
        for w in points.windows(2) {
            run.push(run[run.len() - 1] + w[0].distance(w[1]));
        }
        let total = run[run.len() - 1];
        points
            .into_iter()
            .zip(run.into_iter().map(|l| l / total))
            .collect()
    };
    let fraction = |p: Point| {
        samples
            .iter()
            .min_by(|a, b| a.0.distance(p).total_cmp(&b.0.distance(p)))
            .unwrap()
            .1
    };

    let rounded =
        fillet_edges_variable(&mut model, &joined, std::slice::from_ref(&seam), &hump(), T)
            .unwrap();
    valid(&model, &rounded.shape);
    assert!(rounded.history.is_deleted(&seam));

    let on_upright = |p: &Point| (p.x.hypot(p.y) - 5.0).abs() < 1e-3;
    let on_across = |p: &Point| (p.y.hypot(p.z - 10.0) - 5.0).abs() < 1e-3;
    let near_a_pole = |p: &Point| p.x.hypot(p.z - 10.0) < 1.5;
    let (mut held, mut largest) = (0, 0.0_f64);
    for e in explore_unique(&model, &rounded.shape, ShapeType::Edge).unwrap() {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(&model, &e, T).unwrap() else {
            continue;
        };
        let points: Vec<Point> = (0..=32)
            .map(|i| c.point_at(a + (b - a) * f64::from(i) / 32.0, T).unwrap())
            .collect();
        // A rail lies on one drum, a little off the other along all its run,
        // and bends with the seam: a drum's own seam line is straight.
        let (first, last) = (points[0], points[points.len() - 1]);
        let chord = (last - first) / first.distance(last).max(1e-12);
        if points
            .iter()
            .all(|p| ((*p - first) - chord * (*p - first).dot(chord)).magnitude() < 1e-6)
        {
            continue;
        }
        let upright_rail = points.iter().all(on_upright)
            && points
                .iter()
                .all(|p| (p.y.hypot(p.z - 10.0) - 5.0).abs() < 3.0)
            && points.iter().any(|p| !on_across(p));
        let across_rail = points.iter().all(on_across)
            && points.iter().all(|p| (p.x.hypot(p.y) - 5.0).abs() < 3.0)
            && points.iter().any(|p| !on_upright(p));
        if !upright_rail && !across_rail {
            continue;
        }
        for p in points.iter().filter(|p| !near_a_pole(p)) {
            let (r, other) = if upright_rail {
                let k = (p.z - 10.0).abs() / p.x.abs();
                let centre = Point::new(p.x * k, p.y * k, p.z);
                let s = 5.0 / (5.0 * k);
                let other = Point::new(centre.x, centre.y * s, 10.0 + (centre.z - 10.0) * s);
                (5.0 * k - 5.0, other)
            } else {
                let k = p.x.abs() / (p.z - 10.0).abs();
                let centre = Point::new(p.x, p.y * k, 10.0 + (p.z - 10.0) * k);
                let s = 1.0 / k;
                (
                    5.0 * k - 5.0,
                    Point::new(centre.x * s, centre.y * s, centre.z),
                )
            };
            let off = (r - hump_open(fraction(*p)))
                .abs()
                .min((r - hump_open(fraction(other))).abs());
            assert!(
                off < 1e-3,
                "a rail point {p:?} holds radius {r}, {off} off the law"
            );
            largest = largest.max(r);
            held += 1;
        }
    }
    assert!(held > 32, "{held} rail points");
    assert!(largest > 1.5, "the ball grows to {largest} along the seam");
}

#[test]
fn a_variable_radius_refuses_by_name_what_it_does_not_round() {
    let mut model = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, at(0.0, 0.0, 0.0), 10.0, 20.0, T)
        .unwrap()
        .shape;
    let rim = edge_near(&model, &drum, [0.0, -10.0, 20.0]);
    let linear = RadiusLaw::Linear {
        start: 1.0,
        end: 2.0,
    };
    let refusal = |r: ogeom_core::OgeomResult<ogeom_algo::Built>| match r {
        Err(e) => e.to_string(),
        Ok(_) => panic!("rounded what it should refuse"),
    };
    let jump = refusal(fillet_edges_variable(
        &mut model,
        &drum,
        std::slice::from_ref(&rim),
        &linear,
        T,
    ));
    assert!(
        jump.contains("jumps where a closed chain meets its start"),
        "{jump}"
    );

    let block = ogeom_algo::make_box(&mut model, at(0.0, 0.0, 0.0), (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let two = [
        edge_near(&model, &block, [10.0, 0.0, 10.0]),
        edge_near(&model, &block, [0.0, 10.0, 10.0]),
    ];
    let chain = refusal(fillet_edges_variable(&mut model, &block, &two, &hump(), T));
    assert!(
        chain.contains("open chain of 2 edges is still owed unless"),
        "{chain}"
    );
    for (law, says) in [
        (RadiusLaw::Points(vec![]), "no points"),
        (
            RadiusLaw::Points(vec![(0.5, 1.0), (0.25, 2.0)]),
            "must increase",
        ),
        (RadiusLaw::Points(vec![(1.5, 1.0)]), "off its chain"),
        (
            RadiusLaw::Points(vec![(0.0, 1.0), (1.0, -1.0)]),
            "rounds nothing",
        ),
    ] {
        let why = refusal(fillet_edges_variable(
            &mut model,
            &block,
            &two[..1],
            &law,
            T,
        ));
        assert!(why.contains(says), "{why}");
    }
}
