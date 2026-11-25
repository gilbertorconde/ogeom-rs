//! A convex and a concave round meeting at one vertex: the L-bracket's
//! front edge at its re-entrant edge, (5, 0, 5).
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
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
/// 2, 1 and 1.5, the upright's front edge left sharp (ADR 0001, §13
/// question 6). Both orders give the same solid to the bit, valid, within
/// the bracket and the concave round's fill; patches close the opening
/// between the bands' ends, each patch edge on a face at the vertex or on
/// a band, every patch tangent to its neighbours but along the edge left
/// sharp, which stays sharp up to the vertex.
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
        let past = outside_bracket(&m, &rounded, radius);
        assert!(
            past <= 2e-4,
            "flip {flip}: the corner stands {past} outside the bracket and its fill"
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
        assert!(
            !patches.is_empty(),
            "flip {flip}: no patch closes the corner"
        );
        let mut worst = 0.0_f64;
        for patch in &patches {
            worst = worst.max(tangent_but_at_the_sharp_edge(&m, &rounded, patch));
            // Each edge of the patches' outline lies on a face at the vertex
            // or on a band: the patches close the opening between the
            // bands' ends.
            for e in explore_unique(&m, patch, ShapeType::Edge).unwrap() {
                let between_patches = ogeom_topo::ancestors_of(&m, &rounded, &e, ShapeType::Face)
                    .unwrap()
                    .iter()
                    .all(|f| patches.iter().any(|q| q.is_same(f)));
                if between_patches {
                    continue;
                }
                // A side shrunk to a point has no curve.
                let Some((c, (lo, hi))) = ogeom_algo::placed_edge_curve(&m, &e, T).unwrap() else {
                    continue;
                };
                let p = c.point_at(f64::midpoint(lo, hi), T).unwrap();
                assert!(
                    on_face_or_band(p, radius),
                    "flip {flip}: a patch edge at {p:?}"
                );
            }
        }
        eprintln!(
            "flip {flip}: the patches meet their neighbours within {worst} rad, volume {}",
            volume(&m, &rounded)
        );
        assert!(worst < 2e-4, "flip {flip}: tangent within {worst} rad");
    }
}

/// All three edges at the vertex rounded at 2: one patch closes the
/// opening between the bands' ends, tangent to every band it meets, each of
/// its edges on a face at the vertex or on a band. Every order of
/// the request gives the same solid, to the bit. Outside a cell round the
/// corner, each band is the straight prism of its section, and the solid
/// there measures exactly what the bracket and the three sections give.
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
        let past = outside_bracket(&m, &rounded, 2.0);
        assert!(
            past <= 2e-4,
            "{order:?}: the corner stands {past} outside the bracket and its fill"
        );
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
        assert!(
            contacts.len() >= 3,
            "the patch meets at least the three band ends: {contacts:?}"
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
            // A side shrunk to a point has no curve.
            let Some((c, (lo, hi))) = ogeom_algo::placed_edge_curve(&m, &e, T).unwrap() else {
                continue;
            };
            let p = c.point_at(f64::midpoint(lo, hi), T).unwrap();
            assert!(on_face_or_band(p, 2.0), "a patch edge at {p:?}");
        }
        // Outside the cell, each band is the straight prism of its section,
        // so the solid there is exact: the bracket outside the cell, less the
        // front band and the upright's band and plus the re-entrant band's
        // fill, each section `r^2 (1 - pi / 4)` times its length outside the
        // cell (the front band from x = 9.25 to 20, the fill from y = 4.25 to
        // 20, the upright's band from z = 9.25 to 20).
        // The cell, a quarter wider all round, so its walls stand off the
        // bands' section planes the patch meets them in.
        let cell = ogeom_algo::make_box(
            &mut m,
            Frame::new(Point::new(0.75, -4.25, 0.75), Direction::Z, Direction::X, T).unwrap(),
            (8.5, 8.5, 8.5),
            T,
        )
        .unwrap()
        .shape;
        let left = ogeom_bool::cut(&mut m, &rounded, &cell, T).unwrap().shape;
        let got = volume(&m, &left);
        let section = 4.0 - core::f64::consts::PI;
        // The bracket within the cell: the base there, and the upright above
        // the base's top.
        let within = 4.25 * 4.25 * (8.5 + 4.25);
        let want = (3500.0 - within) - section * (20.0 - 9.25) + section * (20.0 - 4.25)
            - section * (20.0 - 9.25);
        assert!(
            (got - want).abs() < 1e-9 * want,
            "outside the cell {got} against {want}"
        );
    }
}

/// Bands of radius 3 on an upright 5 thick, with all three edges rounded or
/// the upright's front edge left sharp. The rounds stay clear of the
/// upright's far face, so the corner closes as it does at
/// smaller radii (ADR 0001, §13 question 7: no corner of mixed convexity is
/// refused). In either order the solid is the same to the bit and valid,
/// it lies within the bracket and the concave round's fill, and every
/// fitted face is tangent to its neighbours.
#[test]
fn a_signed_corner_as_wide_as_the_wall_closes_in_either_order() {
    for three in [false, true] {
        let mut first: Option<(Vec<[u64; 3]>, usize, u64)> = None;
        for flip in [false, true] {
            let mut m = Model::new();
            let l = bracket(&mut m);
            let [front, crease, upright] = corner_edges(&m, &l);
            let mut edges = vec![front, crease];
            if three {
                edges.push(upright);
            }
            if flip {
                edges.reverse();
            }
            let what = format!("three {three}, flip {flip}");
            let rounded = ogeom_fillet::fillet_edges(&mut m, &l, &edges, 3.0, T)
                .unwrap_or_else(|e| panic!("{what}: {e}"))
                .shape;
            let diagnosis = ogeom_algo::check(&m, &rounded, T).unwrap();
            assert!(diagnosis.is_valid(), "{what}: {:?}", diagnosis.problems);
            let past = outside_bracket(&m, &rounded, 3.0);
            assert!(
                past <= 2e-4,
                "{what}: the corner stands {past} outside the bracket and its fill"
            );
            let faces = explore_unique(&m, &rounded, ShapeType::Face).unwrap();
            let mut worst = 0.0_f64;
            for face in &faces {
                let d = m.node(face).unwrap().data().as_face().unwrap();
                if !matches!(
                    m.geometry().surface(d.surface),
                    Some(ogeom_geom::SurfaceGeometry::BSpline(_))
                ) {
                    continue;
                }
                worst = worst.max(tangent_but_at_the_sharp_edge(&m, &rounded, face));
            }
            assert!(worst < 1e-3, "{what}: tangent within {worst} rad");
            let got = (
                corners_of(&m, &rounded),
                faces.len(),
                volume(&m, &rounded).to_bits(),
            );
            match &first {
                None => first = Some(got),
                Some(want) => assert!(*want == got, "{what} lands on another solid"),
            }
        }
    }
}

/// The largest angle `patch` meets a neighbour at, sampled along each edge
/// they share, but along the upright's front edge's line (x = 5, y = 0):
/// an edge left sharp stays sharp up to the vertex, between whatever lies
/// either side of it there.
fn tangent_but_at_the_sharp_edge(m: &Model, shape: &Shape, patch: &Shape) -> f64 {
    let mut worst = 0.0_f64;
    for c in ogeom_fillet::analyse_blend(m, shape, patch, 17, T).unwrap() {
        let on_sharp = ogeom_algo::placed_edge_curve(m, &c.edge, T)
            .unwrap()
            .is_some_and(|(curve, (lo, hi))| {
                [lo, f64::midpoint(lo, hi), hi].iter().all(|t| {
                    let p = curve.point_at(*t, T).unwrap();
                    (p.x - 5.0).abs() < 1e-6 && p.y.abs() < 1e-6
                })
            });
        if !on_sharp {
            worst = worst.max(c.tangency_error);
        }
    }
    worst
}

/// Whether `p` lies on a face of the bracket at the vertex (x = 5, y = 0
/// or z = 5) or on the surface of one of the three bands of radius `r`
/// there: the front edge's cylinder about y = r, z = 5 - r, the re-entrant
/// edge's about x = z = 5 + r, the upright's front edge's about x = 5 - r,
/// y = r.
fn on_face_or_band(p: Point, r: f64) -> bool {
    let near = |d: f64| d.abs() < 1e-6;
    near(p.x - 5.0)
        || near(p.y)
        || near(p.z - 5.0)
        || near((p.y - r).hypot(p.z - (5.0 - r)) - r)
        || near((p.x - (5.0 + r)).hypot(p.z - (5.0 + r)) - r)
        || near((p.x - (5.0 - r)).hypot(p.y - r) - r)
}

/// How far any point of `shape`'s mesh within 10 of the vertex (5, 0, 5)
/// stands outside the bracket together with the concave round of radius
/// `r` on its re-entrant edge (the fill between the base's top, the
/// upright's face and the round's ball). A convex round only takes
/// material off, so a rounded corner lies within that.
fn outside_bracket(m: &Model, shape: &Shape, r: f64) -> f64 {
    let deflection = ogeom_mesh::Deflection {
        chord: 1e-3,
        ..ogeom_mesh::Deflection::default()
    };
    let out = |p: Point| -> f64 {
        let in_box = |lo: [f64; 3], hi: [f64; 3]| {
            let q = [p.x, p.y, p.z];
            (0..3)
                .map(|i| (lo[i] - q[i]).max(q[i] - hi[i]))
                .fold(f64::NEG_INFINITY, f64::max)
        };
        let base = in_box([0.0, 0.0, 0.0], [20.0, 20.0, 5.0]);
        let upright = in_box([0.0, 0.0, 0.0], [5.0, 20.0, 20.0]);
        let fill = in_box([5.0, 0.0, 5.0], [5.0 + r, 20.0, 5.0 + r])
            .max(r - (p.x - 5.0 - r).hypot(p.z - 5.0 - r));
        base.min(upright).min(fill)
    };
    let mut worst = f64::NEG_INFINITY;
    for f in explore_unique(m, shape, ShapeType::Face).unwrap() {
        for p in ogeom_mesh::triangulate_face(m, &f, deflection, T)
            .unwrap()
            .positions
        {
            if p.distance(Point::new(5.0, 0.0, 5.0)) <= 10.0 {
                worst = worst.max(out(p));
            }
        }
    }
    worst
}

/// What the bracket and its concave round hold near (5, 0, 5), as the
/// region the engine proves a corner patch within: the front face's half
/// space, met with the base's top joined to the upright's front (the
/// re-entrant edge is concave), joined with the concave round's fill there
/// between its two faces, off its ball (the cylinder of radius `r` about
/// (5 + r, y, 5 + r)), short of the chord between the ball's contacts and
/// along the edge from the vertex. A convex round only takes material off.
/// Every piece's height is raised by `shift`: the region shrunk by it.
fn held_near_corner(r: f64, shift: f64) -> ogeom_algo::HalfSpaces {
    use ogeom_algo::HalfSpaces::{Any, Below, Every, Revolved};
    use ogeom_math::Vector;
    let corner = Point::new(5.0, 0.0, 5.0);
    let below = |normal: Vector, at: Point| Below {
        normal,
        at: at - normal * shift,
    };
    let (x, y, z) = (Vector::X, Vector::Y, Vector::Z);
    let solid = Every(vec![
        below(-y, corner),
        Any(vec![below(z, corner), below(x, corner)]),
    ]);
    let towards = (x + z) * (1.0 / 2.0_f64.sqrt());
    let fill = Every(vec![
        below(-z, corner),
        below(-x, corner),
        Revolved {
            at: Point::new(5.0 + r, 0.0, 5.0 + r),
            axis: y,
            radius: r + shift,
            slope: 0.0,
            inside: false,
        },
        below(towards, Point::new(5.0 + r, 0.0, 5.0)),
        below(-y, corner),
    ]);
    Any(vec![solid, fill])
}

/// The engine's corner patch at the bracket's signed corner, with and
/// without the upright's edge rounded, proved within what the solid and its
/// concave round hold, the proof cross-checked by a fine mesh: both clear
/// every patch at the bands' fit target; asked to keep the patch a tenth of
/// that target inside the region (the region shrunk by the target and a
/// tenth, the proof then held to the target), which its sides on the faces
/// cannot, the mesh finds a point past it and the proof clears nothing the
/// mesh does not.
#[test]
fn the_signed_corner_s_patch_is_proved_within_what_the_solid_holds_as_a_fine_mesh_finds() {
    let fit = 2e-4;
    let r = 2.0;
    for three in [false, true] {
        let mut m = Model::new();
        let l = bracket(&mut m);
        let [front, inner, upright] = corner_edges(&m, &l);
        let mut edges = vec![front, inner];
        if three {
            edges.push(upright);
        }
        let rounded = ogeom_fillet::fillet_edges(&mut m, &l, &edges, r, T)
            .unwrap()
            .shape;
        let patches: Vec<Shape> = explore_unique(&m, &rounded, ShapeType::Face)
            .unwrap()
            .into_iter()
            .filter(|f| {
                let data = m.node(f).unwrap().data().as_face().unwrap();
                matches!(
                    m.geometry().surface(data.surface),
                    Some(ogeom_geom::SurfaceGeometry::BSpline(_))
                ) && ogeom_algo::shape_bounds(&m, f, T)
                    .unwrap()
                    .distance_to(Point::new(5.0, 0.0, 5.0))
                    < 2.0 * r
            })
            .collect();
        assert!(!patches.is_empty(), "three {three}: the corner's patches");
        let region = held_near_corner(r, 0.0);
        let shrunk = held_near_corner(r, 1.1 * fit);
        let fine = ogeom_mesh::Deflection {
            chord: fit * 0.25,
            ..ogeom_mesh::Deflection::default()
        };
        for patch in &patches {
            let mesh = ogeom_mesh::triangulate_face(&m, patch, fine, T).unwrap();
            let sampled = mesh
                .positions
                .iter()
                .map(|p| region.height(*p).0)
                .fold(f64::NEG_INFINITY, f64::max);
            let proved = ogeom_algo::certify_face_within(&m, patch, &region, fit, T).unwrap();
            assert!(
                sampled <= fit,
                "three {three}: the mesh finds {sampled} outside"
            );
            assert!(
                proved.is_ok(),
                "three {three}: the proof refuses {proved:?}"
            );
            let strict = -0.1 * fit;
            let proved = ogeom_algo::certify_face_within(&m, patch, &shrunk, fit, T).unwrap();
            assert!(
                sampled > strict,
                "three {three}: the mesh finds every point well inside"
            );
            assert!(
                proved.is_err(),
                "three {three}: the proof clears what the mesh does not"
            );
        }
    }
}
