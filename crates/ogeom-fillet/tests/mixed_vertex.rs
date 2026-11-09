//! A convex and a concave round meeting at one vertex: the L-bracket's
//! front edge at its re-entrant edge, (5, 0, 5).
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn edge_near(m: &Model, s: &Shape, q: Point) -> Shape {
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

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// The L-bracket: a base 20 x 20 x 5 and an upright 5 x 20 x 20.
fn bracket(m: &mut Model) -> Shape {
    let bx = |m: &mut Model, s: (f64, f64, f64)| {
        ogeom_algo::make_box(
            m,
            Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap(),
            s,
            T,
        )
        .unwrap()
        .shape
    };
    let a = bx(m, (20.0, 20.0, 5.0));
    let b = bx(m, (5.0, 20.0, 20.0));
    ogeom_bool::fuse(m, &a, &b, T).unwrap().shape
}

/// The bracket's three edges at (5, 0, 5): the front edge along the base's
/// top (convex), the re-entrant edge (concave) and the upright's front
/// edge (convex).
fn corner_edges(m: &Model, l: &Shape) -> [Shape; 3] {
    [
        edge_near(m, l, Point::new(12.5, 0.0, 5.0)),
        edge_near(m, l, Point::new(5.0, 10.0, 5.0)),
        edge_near(m, l, Point::new(5.0, 0.0, 12.5)),
    ]
}

/// Every vertex of `s`, in one order, to the bit.
fn corners_of(m: &Model, s: &Shape) -> Vec<[u64; 3]> {
    let mut all: Vec<[u64; 3]> = explore_unique(m, s, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| {
            let p = m.node(v).unwrap().data().as_vertex().unwrap().point;
            let p = v.transform(m.datums()).unwrap().apply(p);
            [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
        })
        .collect();
    all.sort_unstable();
    all
}

/// With the third edge at the vertex left sharp, the opening between the
/// two bands' ends would cross it, and no one patch is tangent to both of
/// its faces there: the request is refused by name in either order, where
/// it once took whichever corner was asked first.
#[test]
fn a_convex_and_a_concave_round_beside_an_edge_left_sharp_are_refused_in_either_order() {
    for flip in [false, true] {
        let mut m = Model::new();
        let l = bracket(&mut m);
        let [front, crease, _] = corner_edges(&m, &l);
        let edges = if flip {
            vec![crease, front]
        } else {
            vec![front, crease]
        };
        let Err(e) = ogeom_fillet::fillet_edges(&mut m, &l, &edges, 2.0, T) else {
            panic!("flip {flip}: a corner taken in the order asked");
        };
        assert!(
            e.to_string()
                .contains("the request rounds every edge of the vertex"),
            "flip {flip}: {e}"
        );
    }
}

/// All three edges at the vertex rounded at 2: each band (2 wide on its
/// faces) stops 4 from the vertex, on x = 9, y = 4 and z = 9, and one patch
/// closes the opening, tangent to every band end and face. Every order of
/// the request gives the same solid, to the bit. Outside the cell the
/// corner rebuilds (the box from each section back to 4 behind the vertex)
/// the solid is the three bands built one after another.
#[test]
fn three_rounds_of_both_kinds_at_a_vertex_close_one_corner_whatever_the_order() {
    let orders = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let mut first: Option<(Vec<[u64; 3]>, usize, u64)> = None;
    for order in orders {
        let mut m = Model::new();
        let l = bracket(&mut m);
        let edges = corner_edges(&m, &l);
        let asked: Vec<Shape> = order.iter().map(|&k| edges[k].clone()).collect();
        let rounded = ogeom_fillet::fillet_edges(&mut m, &l, &asked, 2.0, T)
            .unwrap_or_else(|e| panic!("{order:?}: {e}"))
            .shape;
        let diagnosis = ogeom_algo::check(&m, &rounded, T).unwrap();
        assert!(diagnosis.is_valid(), "{order:?}: {:?}", diagnosis.problems);
        let faces = explore_unique(&m, &rounded, ShapeType::Face).unwrap().len();
        let got = (
            corners_of(&m, &rounded),
            faces,
            volume(&m, &rounded).to_bits(),
        );
        match &first {
            None => first = Some(got),
            Some(want) => assert!(*want == got, "{order:?} lands on another solid"),
        }
        if order != [0, 1, 2] {
            continue;
        }
        // The patch: the one fitted face.
        let patches: Vec<Shape> = explore_unique(&m, &rounded, ShapeType::Face)
            .unwrap()
            .into_iter()
            .filter(|f| {
                let d = m.node(f).unwrap().data().as_face().unwrap();
                matches!(
                    m.geometry().surface(d.surface),
                    Some(ogeom_geom::SurfaceGeometry::BSpline(_))
                )
            })
            .collect();
        let [patch] = patches.as_slice() else {
            panic!("{} fitted faces", patches.len());
        };
        let contacts = ogeom_fillet::analyse_blend(&m, &rounded, patch, 17, T).unwrap();
        assert_eq!(
            contacts.len(),
            6,
            "three band ends and three faces: {contacts:?}"
        );
        let tangency = contacts
            .iter()
            .map(|c| c.tangency_error)
            .fold(0.0, f64::max);
        eprintln!("the signed corner's patch meets its neighbours within {tangency} rad");
        assert!(
            tangency < 1e-3,
            "the patch meets its neighbours within {tangency} rad"
        );
        for e in explore_unique(&m, patch, ShapeType::Edge).unwrap() {
            let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(&m, &e, T).unwrap().unwrap();
            let p = c.point_at(f64::midpoint(lo, hi), T).unwrap();
            let on_face = (p.x - 5.0).abs() < 1e-6 || p.y.abs() < 1e-6 || (p.z - 5.0).abs() < 1e-6;
            let on_end =
                (p.x - 9.0).abs() < 1e-6 || (p.y - 4.0).abs() < 1e-6 || (p.z - 9.0).abs() < 1e-6;
            assert!(on_face || on_end, "a patch edge at {p:?}");
        }
        // Outside the cell, the bands as one after another builds them.
        let mut stepwise = l.clone();
        for k in order {
            let edge = {
                let [a, b, c] = corner_edges(&m, &stepwise);
                [a, b, c][k].clone()
            };
            stepwise = ogeom_fillet::fillet_edge(&mut m, &stepwise, &edge, 2.0, T)
                .unwrap()
                .shape;
        }
        let cell = ogeom_algo::make_box(
            &mut m,
            Frame::new(Point::new(1.0, -4.0, 1.0), Direction::Z, Direction::X, T).unwrap(),
            (8.0, 8.0, 8.0),
            T,
        )
        .unwrap()
        .shape;
        let outside = |m: &mut Model, s: &Shape| {
            let left = ogeom_bool::cut(m, s, &cell, T).unwrap().shape;
            volume(m, &left)
        };
        let (a, b) = (outside(&mut m, &rounded), outside(&mut m, &stepwise));
        assert!((a - b).abs() < 1e-9, "outside the cell {a} against {b}");
        eprintln!(
            "the signed corner {} against the stepwise bands {}",
            volume(&m, &rounded),
            volume(&m, &stepwise)
        );
    }
}
