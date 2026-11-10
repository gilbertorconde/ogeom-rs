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

/// The issue's request: the front edge and the re-entrant edge rounded at
/// 2, the upright's front edge left sharp. The opening between the two
/// bands' ends crosses the sharp edge; it runs into a point there, where a
/// patch with a shrunk side turns its normal from the front face's to the
/// upright's, and one patch closes the rest. Both orders give the same
/// solid to the bit, valid, every patch tangent to its neighbours, each
/// band stopped short of the vertex on its section plane.
#[test]
fn a_convex_and_a_concave_round_beside_an_edge_left_sharp_close_one_corner_in_either_order() {
    for radius in [2.0, 1.0, 1.5] {
        two_edges_in_either_order(radius);
    }
}

/// The issue's two edges rounded at `radius`, in both orders.
fn two_edges_in_either_order(radius: f64) {
    let mut first: Option<(Vec<[u64; 3]>, usize, u64)> = None;
    for flip in [false, true] {
        let mut m = Model::new();
        let l = bracket(&mut m);
        let [front, crease, _] = corner_edges(&m, &l);
        let edges = if flip {
            vec![crease, front]
        } else {
            vec![front, crease]
        };
        let rounded = ogeom_fillet::fillet_edges(&mut m, &l, &edges, radius, T)
            .unwrap_or_else(|e| panic!("flip {flip}: {e}"))
            .shape;
        let diagnosis = ogeom_algo::check(&m, &rounded, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "flip {flip}: {:?}",
            diagnosis.problems
        );
        let faces = explore_unique(&m, &rounded, ShapeType::Face).unwrap();
        let got = (
            corners_of(&m, &rounded),
            faces.len(),
            volume(&m, &rounded).to_bits(),
        );
        match &first {
            None => first = Some(got),
            Some(want) => assert!(*want == got, "flip {flip} lands on another solid"),
        }
        // The two fitted faces: the patch at the sharp edge and the rest.
        let patches: Vec<Shape> = faces
            .into_iter()
            .filter(|f| {
                let d = m.node(f).unwrap().data().as_face().unwrap();
                matches!(
                    m.geometry().surface(d.surface),
                    Some(ogeom_geom::SurfaceGeometry::BSpline(_))
                )
            })
            .collect();
        assert_eq!(
            patches.len(),
            2,
            "flip {flip}: a patch at the sharp edge and one over the rest"
        );
        let mut worst = 0.0_f64;
        let mut ends = (false, false);
        for patch in &patches {
            for c in ogeom_fillet::analyse_blend(&m, &rounded, patch, 17, T).unwrap() {
                worst = worst.max(c.tangency_error);
            }
            // Each band's end arc lies on its section, short of the vertex:
            // the front band's on a plane x = c past 5, the re-entrant
            // band's on a plane y = c past 0.
            for e in explore_unique(&m, patch, ShapeType::Edge).unwrap() {
                // The patch at the sharp edge's point has a side shrunk to
                // it, an edge with no curve.
                let Some((c, (lo, hi))) = ogeom_algo::placed_edge_curve(&m, &e, T).unwrap() else {
                    continue;
                };
                let p: Vec<Point> = [lo, f64::midpoint(lo, hi), hi]
                    .iter()
                    .map(|t| c.point_at(*t, T).unwrap())
                    .collect();
                let flat = |f: fn(&Point) -> f64| {
                    (f(&p[0]) - f(&p[1])).abs() < 1e-9 && (f(&p[0]) - f(&p[2])).abs() < 1e-9
                };
                if flat(|q| q.x) && p[0].x > 5.0 + 1e-6 && p[1].y > 1e-6 {
                    ends.0 = true;
                }
                if flat(|q| q.y) && p[0].y > 1e-6 && p[1].x > 5.0 + 1e-6 {
                    ends.1 = true;
                }
            }
        }
        eprintln!(
            "flip {flip}: the patches meet their neighbours within {worst} rad, volume {}",
            volume(&m, &rounded)
        );
        assert!(worst < 2e-4, "flip {flip}: tangent within {worst} rad");
        assert!(
            ends.0 && ends.1,
            "flip {flip}: a band ends on the other's rail: {ends:?}"
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

/// Bands as wide as the upright is thick: the cell the corner cuts reaches
/// twice the widest band behind the vertex, through the upright's far
/// face, which nothing would put back. Refused by name, with all three
/// edges rounded or the upright's front edge left sharp, as the cell
/// guard refuses any other feature inside the cell.
#[test]
fn a_signed_corner_whose_cell_reaches_through_the_wall_is_refused_by_name() {
    for three in [false, true] {
        let mut m = Model::new();
        let l = bracket(&mut m);
        let [front, crease, upright] = corner_edges(&m, &l);
        let mut edges = vec![front, crease];
        if three {
            edges.push(upright);
        }
        let err = ogeom_fillet::fillet_edges(&mut m, &l, &edges, 3.0, T)
            .map(|_| ())
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("another feature of the solid lies inside the setback corner's cell"),
            "three {three}: {err}"
        );
    }
}
