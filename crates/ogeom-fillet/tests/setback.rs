//! Setback corners: the bands at a vertex stopped a chosen distance along
//! their edges, and one patch over the opening between them.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{
    Setbacks, analyse_blend, fillet_edges, fillet_edges_cornered, round_vertex_setback,
};
use ogeom_geom::{Continuity, Curve3d as _};
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

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

/// The vertex of `s` standing at `p`.
fn vertex_at(m: &Model, s: &Shape, p: [f64; 3]) -> Shape {
    let q = Point::new(p[0], p[1], p[2]);
    explore_unique(m, s, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .find(|v| {
            m.node(v)
                .and_then(|n| n.data().as_vertex().map(|d| d.point))
                .is_some_and(|x| x.distance(q) < 1e-9)
        })
        .expect("a vertex there")
}

fn volume(m: &Model, s: &Shape) -> f64 {
    let deflection = ogeom_mesh::Deflection {
        chord: 1e-3,
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

/// A 20 x 20 x 10 block and its corner at the origin's top, (0, 0, 10),
/// with the three edges leaving it.
fn block() -> (Model, Shape, Shape, Vec<Shape>) {
    let mut m = Model::new();
    let b = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges: Vec<Shape> = [[10.0, 0.0, 10.0], [0.0, 10.0, 10.0], [0.0, 0.0, 5.0]]
        .iter()
        .map(|p| edge_near(&m, &b, *p))
        .collect();
    let v = vertex_at(&m, &b, [0.0, 0.0, 10.0]);
    (m, b, v, edges)
}

fn even(setback: f64, continuity: Continuity) -> Setbacks {
    Setbacks {
        per_edge: Vec::new(),
        default: setback,
        continuity,
    }
}

#[test]
fn a_setback_of_the_radius_is_the_ball_s_corner() {
    let (mut m, b, v, edges) = block();
    let setback = round_vertex_setback(&mut m, &b, &v, 2.0, &even(2.0, Continuity::G1), T)
        .unwrap()
        .shape;
    let ball = fillet_edges(&mut m, &b, &edges, 2.0, T).unwrap().shape;
    valid(&m, &setback);
    let (a, c) = (volume(&m, &setback), volume(&m, &ball));
    assert!((a - c).abs() < 1e-12, "{a} against {c}");
}

/// The patch a setback corner of `continuity` builds at the block's
/// corner, its bands stopped 5 from the vertex, and the patch's worst
/// tangency against the faces round it.
fn five_back(continuity: Continuity) -> (Model, Shape, Shape, f64) {
    let (mut m, b, v, _) = block();
    let rounded = round_vertex_setback(&mut m, &b, &v, 2.0, &even(5.0, continuity), T).unwrap();
    valid(&m, &rounded.shape);
    let patch = rounded.history.generated(&v)[0].clone();
    let contacts = analyse_blend(&m, &rounded.shape, &patch, 17, T).unwrap();
    assert_eq!(
        contacts.len(),
        6,
        "three band ends and three faces: {contacts:?}"
    );
    let tangency = contacts
        .iter()
        .map(|c| c.tangency_error)
        .fold(0.0, f64::max);
    (m, rounded.shape, patch, tangency)
}

#[test]
fn a_setback_corner_stops_each_band_at_its_setback_tangent_to_the_patch() {
    let (m, rounded, patch, tangency) = five_back(Continuity::G1);
    assert!(
        tangency < 1e-3,
        "the patch meets its neighbours within {tangency} rad"
    );
    // Every point of the patch's boundary on a band's end lies 5 along
    // that band's edge from the vertex: on x = 5, y = 5 or z = 10 - 5.
    for e in explore_unique(&m, &patch, ShapeType::Edge).unwrap() {
        let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(&m, &e, T).unwrap().unwrap();
        let p = c.point_at(f64::midpoint(lo, hi), T).unwrap();
        let on_face = p.x.abs() < 1e-6 || p.y.abs() < 1e-6 || (p.z - 10.0).abs() < 1e-6;
        let on_end =
            (p.x - 5.0).abs() < 1e-6 || (p.y - 5.0).abs() < 1e-6 || (p.z - 5.0).abs() < 1e-6;
        assert!(on_face || on_end, "a patch edge at {p:?}");
    }
    // The setback corner takes off more than the ball's corner does.
    let (mut m2, b, v, edges) = block();
    let ball = fillet_edges(&mut m2, &b, &edges, 2.0, T).unwrap().shape;
    let _ = v;
    assert!(volume(&m, &rounded) < volume(&m2, &ball));
}

#[test]
fn a_setback_corner_asked_g2_is_refused_by_name_while_bands_are_circular() {
    // At a band's arc end on a face the arc bends the patch at 1/r along
    // itself; the face asks no bend across its trim, which runs along
    // the arc's direction there.
    let (mut m, b, v, _) = block();
    let refused = round_vertex_setback(&mut m, &b, &v, 2.0, &even(5.0, Continuity::G2), T);
    let Err(e) = refused else {
        panic!("a G2 setback corner built against circular bands");
    };
    assert!(
        e.to_string()
            .contains("curvature continuous with their faces"),
        "{e}"
    );
}

#[test]
fn a_setback_corner_rounds_a_square_pyramid_s_apex() {
    let mut m = Model::new();
    let base = ogeom_algo::make_box(&mut m, at(-10.0, -10.0, -10.0), (20.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let mut solid = base;
    // Four planes leaning in from the base's top edges to the apex at
    // (0, 0, 10), each cutting a half-space off.
    for (nx, ny) in [(1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)] {
        let normal = ogeom_math::Vector::new(nx, ny, 1.0);
        let across = ogeom_math::Vector::new(-ny, nx, 0.0);
        let side = normal.cross(across) / normal.cross(across).magnitude();
        let corner = Point::new(0.0, 0.0, 10.0) - across * 40.0 - side * 40.0;
        let frame = Frame::new(
            corner,
            Direction::new(normal, T).unwrap(),
            Direction::new(across, T).unwrap(),
            T,
        )
        .unwrap();
        let away = ogeom_algo::make_box(&mut m, frame, (80.0, 80.0, 40.0), T)
            .unwrap()
            .shape;
        solid = ogeom_bool::cut(&mut m, &solid, &away, T).unwrap().shape;
    }
    let apex = vertex_at(&m, &solid, [0.0, 0.0, 10.0]);
    let rounded =
        round_vertex_setback(&mut m, &solid, &apex, 1.0, &even(4.0, Continuity::G1), T).unwrap();
    valid(&m, &rounded.shape);
}

#[test]
fn a_setback_shorter_than_the_ball_s_is_refused() {
    let (mut m, b, v, _) = block();
    let refused = round_vertex_setback(&mut m, &b, &v, 2.0, &even(1.5, Continuity::G1), T);
    let Err(e) = refused else {
        panic!("a setback inside the ball's corner rounded");
    };
    assert!(
        e.to_string().contains("inside the ball's own corner"),
        "{e}"
    );
}

/// A sharp corner's case: the edge left sharp, the radius, the setbacks
/// named by edge index, and the default.
type SharpCase = (usize, f64, &'static [(usize, f64)], f64);

/// The block's corner with edge `sharp` of [`block`]'s left sharp and the
/// other two rounded at `radius`, set back as `per_edge` names by index
/// and `default` otherwise: valid, two patches, the
/// worst tangency of either against its neighbours, and, along the shrunk
/// side of the patch at the sharp edge, how its normal turns.
fn sharp_corner(
    sharp: usize,
    radius: f64,
    per_edge: &[(usize, f64)],
    default: f64,
) -> (f64, Vec<f64>, f64) {
    use ogeom_geom::Surface as _;
    let (mut m, b, v, edges) = block();
    let setbacks = Setbacks {
        per_edge: per_edge
            .iter()
            .map(|(k, s)| (edges[*k].clone(), *s))
            .collect(),
        default,
        continuity: Continuity::G1,
    };
    let rounded = ogeom_fillet::round_vertex_setback_sharp(
        &mut m,
        &b,
        &v,
        &edges[sharp],
        radius,
        &setbacks,
        T,
    )
    .unwrap();
    valid(&m, &rounded.shape);
    // Nothing of it stands outside the block.
    let past = farthest_out(&m, &rounded.shape, |p| {
        [-p.x, -p.y, -p.z, p.x - 20.0, p.y - 20.0, p.z - 10.0]
            .into_iter()
            .fold(f64::NEG_INFINITY, f64::max)
    });
    assert!(past <= 2e-4, "the corner stands {past} outside the block");
    let patches = rounded.history.generated(&v);
    assert_eq!(
        patches.len(),
        2,
        "a patch at the sharp edge and one over the rest"
    );
    let mut worst = 0.0_f64;
    for patch in patches {
        for c in analyse_blend(&m, &rounded.shape, patch, 17, T).unwrap() {
            worst = worst.max(c.tangency_error);
        }
    }
    let data = m.node(&patches[0]).unwrap().data().as_face().unwrap();
    let ogeom_geom::SurfaceGeometry::BSpline(surface) = m.geometry().surface(data.surface).unwrap()
    else {
        panic!("the patch at the sharp edge is a B-spline");
    };
    let turns: Vec<f64> = (0..=32)
        .map(|j| {
            let (du, dv) = surface.d1_at(1e-7, f64::from(j) / 32.0, T).unwrap();
            let n = du.cross(dv);
            let n = n / n.magnitude();
            let (start, _) = surface.d1_at(1e-7, 0.0, T).unwrap();
            let (_, dv0) = surface.d1_at(1e-7, 0.0, T).unwrap();
            let n0 = start.cross(dv0);
            (n.dot(n0 / n0.magnitude())).clamp(-1.0, 1.0).acos()
        })
        .collect();
    (worst, turns, volume(&m, &rounded.shape))
}

/// The block's corner with two edges rounded and the third left sharp: the
/// patch runs into the sharp edge at its setback, a three-sided patch there
/// turning its normal from one of the edge's faces to the other's, and
/// both patches meet every band end and face tangent.
#[test]
fn a_setback_corner_leaving_an_edge_sharp_runs_its_patch_into_that_edge() {
    let cases: [SharpCase; 6] = [
        // The top edges rounded, the upright left sharp, all set back 5.
        (2, 2.0, &[], 5.0),
        // The sharp edge cut farther off.
        (2, 2.0, &[(2, 8.0)], 5.0),
        // Uneven bands: at 6 and 4.
        (2, 2.0, &[(0, 6.0), (1, 4.0)], 5.0),
        // A top edge left sharp, the other top edge and the upright rounded.
        (0, 1.5, &[], 4.0),
        // A smaller and a larger radius.
        (2, 1.0, &[], 3.0),
        (2, 3.0, &[], 8.0),
    ];
    for (sharp, radius, per_edge, default) in cases {
        let (worst, turns, volume) = sharp_corner(sharp, radius, per_edge, default);
        eprintln!(
            "MEASURE sharp {sharp} worst {worst} volume {volume} turns {:?}",
            turns.iter().step_by(8).collect::<Vec<_>>()
        );
        assert!(
            worst < 2e-4,
            "the patches meet their neighbours within {worst} rad"
        );
        // The normal turns from one face's to the other's, a right angle,
        // one way all along the shrunk side.
        let last = turns[turns.len() - 1];
        assert!(
            (last - core::f64::consts::FRAC_PI_2).abs() < 1e-6,
            "turns {last}"
        );
        assert!(turns.windows(2).all(|w| w[1] >= w[0] - 1e-9), "{turns:?}");
        assert!(volume < 4000.0);
    }
}

/// Corners whose patch, tangent to every band end and face and fitted
/// freely, stands past a face of the block (by 0.01 to 0.11): each is
/// either closed with its patches held inside the vertex's planes, valid,
/// tangent to every neighbour and within the block, or refused by name;
/// never built outside the solid.
#[test]
fn a_sharp_corner_whose_free_patch_would_leave_the_block_closes_within_it_or_is_refused() {
    let cases: [SharpCase; 4] = [
        (2, 2.0, &[(2, 3.0)], 5.0),
        (2, 2.0, &[(0, 6.0), (1, 4.0), (2, 3.0)], 5.0),
        (1, 1.0, &[(0, 3.0), (2, 5.0)], 4.0),
        (0, 2.0, &[(1, 7.0)], 5.0),
    ];
    let mut built = 0;
    for (sharp, radius, per_edge, default) in cases {
        let (mut m, b, v, edges) = block();
        let setbacks = Setbacks {
            per_edge: per_edge
                .iter()
                .map(|(k, s)| (edges[*k].clone(), *s))
                .collect(),
            default,
            continuity: Continuity::G1,
        };
        let got = ogeom_fillet::round_vertex_setback_sharp(
            &mut m,
            &b,
            &v,
            &edges[sharp],
            radius,
            &setbacks,
            T,
        );
        match got {
            Ok(rounded) => {
                valid(&m, &rounded.shape);
                let mut worst = 0.0_f64;
                let mut past = f64::NEG_INFINITY;
                for patch in rounded.history.generated(&v) {
                    for c in analyse_blend(&m, &rounded.shape, patch, 17, T).unwrap() {
                        worst = worst.max(c.tangency_error);
                    }
                    let mesh = ogeom_mesh::triangulate_face(
                        &m,
                        patch,
                        ogeom_mesh::Deflection {
                            chord: 2e-5,
                            ..ogeom_mesh::Deflection::default()
                        },
                        T,
                    )
                    .unwrap();
                    for p in &mesh.positions {
                        past = past.max((-p.x).max(-p.y).max(p.z - 10.0));
                    }
                }
                eprintln!(
                    "{sharp} {radius} {per_edge:?}: built, tangent within {worst} rad, {past} past the block"
                );
                assert!(
                    worst < 2e-4,
                    "{sharp} {radius} {per_edge:?}: tangent within {worst}"
                );
                assert!(
                    past <= 2e-4,
                    "{sharp} {radius} {per_edge:?}: {past} outside the block"
                );
                built += 1;
            }
            Err(ogeom_core::OgeomError::NotDone(why)) => {
                eprintln!("{sharp} {radius} {per_edge:?}: refused: {why}");
            }
            Err(e) => panic!("{sharp} {radius} {per_edge:?}: refused otherwise: {e}"),
        }
    }
    eprintln!("{built} of 4 built within the block");
}

/// A sharp edge cut no farther from the vertex than the bands' rails stand
/// from its faces leaves no room for the trims from it: refused by name.
#[test]
fn a_sharp_edge_cut_inside_the_bands_is_refused_by_name() {
    let (mut m, b, v, edges) = block();
    let setbacks = Setbacks {
        per_edge: vec![(edges[2].clone(), 1.5)],
        default: 5.0,
        continuity: Continuity::G1,
    };
    let err =
        ogeom_fillet::round_vertex_setback_sharp(&mut m, &b, &v, &edges[2], 2.0, &setbacks, T)
            .unwrap_err()
            .to_string();
    assert!(err.contains("edge left sharp"), "{err}");
}

/// The block with all twelve edges to round, and its two top corners on
/// the edge y = 0, z = 10: (0, 0, 10) and (20, 0, 10).
fn two_corners() -> (Model, Shape, Vec<Shape>, [Shape; 2]) {
    let (m, b, _, _) = block();
    let edges = explore_unique(&m, &b, ShapeType::Edge).unwrap();
    let corners = [
        vertex_at(&m, &b, [0.0, 0.0, 10.0]),
        vertex_at(&m, &b, [20.0, 0.0, 10.0]),
    ];
    (m, b, edges, corners)
}

#[test]
fn two_setback_corners_on_one_edge_each_take_off_what_it_does_alone() {
    let (mut m, b, edges, [v0, v1]) = two_corners();
    let ball = fillet_edges(&mut m, &b, &edges, 2.0, T).unwrap().shape;
    let cornered = |m: &mut Model, at: &[&Shape]| {
        let corners: Vec<(Shape, Setbacks)> = at
            .iter()
            .map(|v| ((*v).clone(), even(5.0, Continuity::G1)))
            .collect();
        fillet_edges_cornered(m, &b, &edges, 2.0, &corners, T).unwrap()
    };
    let first = cornered(&mut m, &[&v0]).shape;
    let second = cornered(&mut m, &[&v1]).shape;
    let both = cornered(&mut m, &[&v0, &v1]);
    let swapped = cornered(&mut m, &[&v1, &v0]).shape;
    valid(&m, &both.shape);
    valid(&m, &swapped);
    for v in [&v0, &v1] {
        let patch = both.history.generated(v)[0].clone();
        let contacts = analyse_blend(&m, &both.shape, &patch, 17, T).unwrap();
        assert_eq!(contacts.len(), 6, "{contacts:?}");
        let tangency = contacts
            .iter()
            .map(|c| c.tangency_error)
            .fold(0.0, f64::max);
        assert!(
            tangency < 1e-3,
            "a patch meets its neighbours within {tangency} rad"
        );
    }
    // The cells lie 10 apart along the shared edge, so each corner takes
    // off what it takes off alone, in either order.
    let vb = volume(&m, &ball);
    let alone = [vb - volume(&m, &first), vb - volume(&m, &second)];
    assert!(alone.iter().all(|a| *a > 0.0), "{alone:?}");
    for together in [vb - volume(&m, &both.shape), vb - volume(&m, &swapped)] {
        assert!(
            (together - alone[0] - alone[1]).abs() < 1e-9 * vb,
            "two corners take off {together}, each alone {alone:?}"
        );
    }
}

#[test]
fn setback_corners_at_the_radius_are_the_balls_corners() {
    let (mut m, b, edges, [v0, v1]) = two_corners();
    let ball = fillet_edges(&mut m, &b, &edges, 2.0, T).unwrap().shape;
    let cornered = fillet_edges_cornered(
        &mut m,
        &b,
        &edges,
        2.0,
        &[
            (v0, even(2.0, Continuity::G1)),
            (v1, even(2.0, Continuity::G1)),
        ],
        T,
    )
    .unwrap()
    .shape;
    valid(&m, &cornered);
    let (a, c) = (volume(&m, &cornered), volume(&m, &ball));
    assert!((a - c).abs() < 1e-12, "{a} against {c}");
}

#[test]
fn setback_corners_meeting_along_an_edge_or_missing_an_edge_are_refused_by_name() {
    let (mut m, b, edges, _) = two_corners();
    // (0, 0, 10) and (0, 0, 0) share the vertical edge, 10 long.
    let top = vertex_at(&m, &b, [0.0, 0.0, 10.0]);
    let bottom = vertex_at(&m, &b, [0.0, 0.0, 0.0]);
    let err = fillet_edges_cornered(
        &mut m,
        &b,
        &edges,
        2.0,
        &[
            (top.clone(), even(5.0, Continuity::G1)),
            (bottom, even(5.0, Continuity::G1)),
        ],
        T,
    )
    .unwrap_err();
    assert!(err.to_string().contains("meet"), "{err}");
    // A setback corner rounds every edge at its vertex.
    let (mut m2, b2, v2, _) = block();
    let only_two: Vec<Shape> = [[10.0, 0.0, 10.0], [0.0, 10.0, 10.0]]
        .iter()
        .map(|p| edge_near(&m2, &b2, *p))
        .collect();
    let err = fillet_edges_cornered(
        &mut m2,
        &b2,
        &only_two,
        2.0,
        &[(v2.clone(), even(5.0, Continuity::G1))],
        T,
    )
    .unwrap_err();
    assert!(err.to_string().contains("not among the edges"), "{err}");
    let err = fillet_edges_cornered(
        &mut m,
        &b,
        &edges,
        2.0,
        &[
            (top.clone(), even(5.0, Continuity::G1)),
            (top, even(5.0, Continuity::G1)),
        ],
        T,
    )
    .unwrap_err();
    assert!(err.to_string().contains("named twice"), "{err}");
}

/// The farthest any point of `shape`'s mesh stands by `out`, a signed
/// distance outside the original solid.
fn farthest_out(m: &Model, shape: &Shape, out: impl Fn(Point) -> f64) -> f64 {
    let deflection = ogeom_mesh::Deflection {
        chord: 1e-3,
        ..ogeom_mesh::Deflection::default()
    };
    let mut worst = f64::NEG_INFINITY;
    for f in explore_unique(m, shape, ShapeType::Face).unwrap() {
        for p in ogeom_mesh::triangulate_face(m, &f, deflection, T)
            .unwrap()
            .positions
        {
            worst = worst.max(out(p));
        }
    }
    worst
}

/// A prism on a parallelogram of sides 20, its second side at `angle`
/// degrees from x, raised along (`tilt`, 0, 10): the solid, the top of the
/// origin's edge, that vertex's edges (along x, along the slanted side, the
/// raised one), and how far a point stands outside the prism.
#[allow(clippy::type_complexity, reason = "test helper")]
fn prism(
    m: &mut Model,
    angle: f64,
    tilt: f64,
) -> (Shape, Shape, [Shape; 3], Box<dyn Fn(Point) -> f64>) {
    use ogeom_math::Vector;
    let a = angle.to_radians();
    let side = Vector::new(a.cos(), a.sin(), 0.0) * 20.0;
    let x = Vector::new(20.0, 0.0, 0.0);
    let pts = [
        Point::ORIGIN,
        Point::ORIGIN + x,
        Point::ORIGIN + x + side,
        Point::ORIGIN + side,
    ];
    let wire = ogeom_algo::make_polygon(m, &pts, true, T).unwrap().shape;
    let plane: ogeom_geom::SurfaceGeometry = ogeom_geom::PlaneSurface::over(
        ogeom_math::Plane::through(Point::ORIGIN, Direction::Z),
        (-1e3, 1e3),
        (-1e3, 1e3),
    )
    .unwrap()
    .into();
    let face = ogeom_algo::make_face(m, plane, &[wire], T).unwrap().shape;
    let up = Vector::new(tilt, 0.0, 10.0);
    let made = ogeom_algo::make_prism(m, &face, up, T).unwrap().shape;
    let solid = explore_unique(m, &made, ShapeType::Solid)
        .unwrap()
        .into_iter()
        .next()
        .unwrap_or(made);
    let top = Point::ORIGIN + up;
    let vertex = explore_unique(m, &solid, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .find(|v| {
            let p = m.node(v).unwrap().data().as_vertex().unwrap().point;
            v.transform(m.datums()).unwrap().apply(p).distance(top) < 1e-7
        })
        .expect("the prism's corner");
    let edges = [
        edge_near(m, &solid, [top.x + 5.0, top.y, top.z]),
        edge_near(
            m,
            &solid,
            [top.x + side.x / 4.0, top.y + side.y / 4.0, top.z],
        ),
        edge_near(m, &solid, [tilt / 2.0, 0.0, 5.0]),
    ];
    // The prism's six planes, each by a point and its outward normal.
    let inner = Point::ORIGIN + x * 0.5 + side * 0.5 + up * 0.5;
    let mut planes: Vec<(Point, Vector)> = vec![
        (Point::ORIGIN, Vector::new(0.0, 0.0, -1.0)),
        (top, Vector::new(0.0, 0.0, 1.0)),
    ];
    for (on, along) in [
        (Point::ORIGIN, x),
        (Point::ORIGIN + side, x),
        (Point::ORIGIN, side),
        (Point::ORIGIN + x, side),
    ] {
        let n = along.cross(up);
        let mut n = n / n.magnitude();
        if (inner - on).dot(n) > 0.0 {
            n = -n;
        }
        planes.push((on, n));
    }
    let out = move |p: Point| {
        planes
            .iter()
            .map(|(o, n)| (p - *o).dot(*n))
            .fold(f64::NEG_INFINITY, f64::max)
    };
    (solid, vertex, edges, Box::new(out))
}

/// The two corners whose patch stood farthest outside its solid: a setback
/// corner on a prism raised aslant, radius 0.5 set back 1, 4 and 2 (0.153
/// out), and one with its upright left sharp where the walls meet at 60
/// degrees, radius 2 set back 5 (0.096 out). Each is refused by name or
/// built with no point of it outside the solid past the bands' target.
#[test]
fn a_setback_corner_never_stands_outside_its_solid() {
    let check = |name: &str,
                 built: ogeom_core::OgeomResult<ogeom_algo::Built>,
                 m: &Model,
                 out: &dyn Fn(Point) -> f64| {
        match built {
            Ok(b) => {
                valid(m, &b.shape);
                let past = farthest_out(m, &b.shape, out);
                assert!(
                    past <= 2e-4,
                    "{name}: the corner stands {past} outside its solid"
                );
            }
            Err(ogeom_core::OgeomError::NotDone(why)) => eprintln!("{name}: refused: {why}"),
            Err(e) => panic!("{name}: refused otherwise: {e}"),
        }
    };
    let mut m = Model::new();
    let (solid, vertex, edges, out) = prism(&mut m, 90.0, 3.0);
    let setbacks = Setbacks {
        per_edge: edges.iter().cloned().zip([1.0, 4.0, 2.0]).collect(),
        default: 5.0,
        continuity: Continuity::G1,
    };
    let built = round_vertex_setback(&mut m, &solid, &vertex, 0.5, &setbacks, T);
    check("aslant", built, &m, &out);

    let mut m = Model::new();
    let (solid, vertex, edges, out) = prism(&mut m, 60.0, 0.0);
    let built = ogeom_fillet::round_vertex_setback_sharp(
        &mut m,
        &solid,
        &vertex,
        &edges[2],
        2.0,
        &even(5.0, Continuity::G1),
        T,
    );
    check("sixty", built, &m, &out);
}
