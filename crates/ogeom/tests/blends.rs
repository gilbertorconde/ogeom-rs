//! What a blend achieved, measured: blends between faces that share no
//! edge, edges whose envelope has no closed form, and the corner where
//! three of them meet.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The edge of `shape` whose midpoint is nearest `near`.
fn edge_near(model: &Model, shape: &Shape, near: Point) -> Shape {
    use ogeom::geom::Curve3d as _;
    // A degenerate edge (a sphere's pole) has no curve and no midpoint.
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            model
                .node(e)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|d| d.curve3d().is_some())
        })
        .min_by(|a, b| {
            let mid = |e: &Shape| {
                let data = model.node(e).unwrap().data().as_edge().unwrap();
                let ogeom::topo::EdgeRepr::Curve3d { curve, range, .. } = data.curve3d().unwrap()
                else {
                    unreachable!()
                };
                model
                    .geometry()
                    .curve(*curve)
                    .unwrap()
                    .point_at(f64::midpoint(range.0, range.1), T)
                    .unwrap()
                    .distance(near)
            };
            mid(a)
                .partial_cmp(&mid(b))
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .expect("some edge")
}

/// The planar face of `shape` whose plane passes through `on` and whose own
/// vertices bracket it.
fn planar_face_at(model: &Model, shape: &Shape, on: Point) -> Shape {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let ogeom::topo::NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            let Some(ogeom::geom::SurfaceGeometry::Plane(plane)) =
                model.geometry().surface(data.surface)
            else {
                return false;
            };
            if plane.plane().distance_to(on).abs() > 1e-9 {
                return false;
            }
            let mut bound = ogeom::math::Aabb::EMPTY;
            for v in explore_unique(model, f, ShapeType::Vertex).unwrap() {
                bound = bound.with_point(model.node(&v).unwrap().data().as_vertex().unwrap().point);
            }
            bound.expanded(1e-6).contains(on)
        })
        .expect("a planar face there")
}

#[test]
fn a_fillet_reports_its_own_tangency_instead_of_claiming_it() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 12.0), T)
        .unwrap()
        .shape;
    let edge = edge_near(&model, &block, Point::new(20.0, 0.0, 12.0));
    let blended = ogeom::fillet::fillet_edge(&mut model, &block, &edge, 2.0, T)
        .unwrap()
        .shape;

    // The blend is the one cylindrical face on the result.
    let blend = explore_unique(&model, &blended, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let ogeom::topo::NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            matches!(
                model.geometry().surface(data.surface),
                Some(ogeom::geom::SurfaceGeometry::Cylinder(_))
            )
        })
        .expect("the rolling ball left a cylinder");

    let contacts = ogeom::fillet::analyse_blend(&model, &blended, &blend, 9, T).unwrap();
    assert_eq!(contacts.len(), 4, "two tangency edges and two end caps");
    // The two long edges are the tangency lines: smooth to rounding. The
    // two ends are the cap arcs, where the blend meets a face it is *not*
    // tangent to: a right angle, and it should say so.
    let mut smooth = 0;
    let mut square = 0;
    for contact in &contacts {
        assert!(
            contact.gap < 1e-9,
            "the shared edge lies on both surfaces: {}",
            contact.gap
        );
        if contact.tangency_error < 1e-9 {
            smooth += 1;
        } else if (contact.tangency_error - core::f64::consts::FRAC_PI_2).abs() < 1e-9 {
            square += 1;
        }
    }
    assert_eq!(
        (smooth, square),
        (2, 2),
        "two tangent joins, two square ones: {contacts:?}"
    );
}

#[test]
fn a_blend_bridges_two_faces_that_share_no_edge() {
    // A step: a tall block and a low one side by side, their vertical wall
    // and horizontal lid meeting at no edge at all. The rolling ball still
    // has a seat (it touches both), and the blend is the fillet that seat
    // implies.
    let mut model = Model::new();
    let tall = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let low = ogeom::algo::make_box(
        &mut model,
        Frame::new(Point::new(10.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
        (20.0, 20.0, 10.0),
        T,
    )
    .unwrap()
    .shape;
    let step = ogeom::boolean::fuse(&mut model, &tall, &low, T)
        .unwrap()
        .shape;

    let wall = planar_face_at(&model, &step, Point::new(10.0, 10.0, 15.0));
    let lid = planar_face_at(&model, &step, Point::new(20.0, 10.0, 10.0));

    let blended = ogeom::fillet::blend_faces(&mut model, &step, &wall, &lid, 4.0, T)
        .unwrap()
        .shape;
    let volume =
        ogeom::algo::volume_properties(&model, &blended, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    // The step is 10*20*20 + 20*20*10 = 8000, and its inner corner is
    // concave: the ball rolls in the notch, so the blend *fills* it with
    // what a square corner would have held minus the quarter disc,
    // (r^2 - pi r^2 / 4), along the 20 of run.
    let r: f64 = 4.0;
    let filled = r.mul_add(r, -(core::f64::consts::PI * r * r / 4.0)) * 20.0;
    assert!(
        (volume - (8000.0 + filled)).abs() < 8000.0 * 2e-3,
        "the notch is filled, not cut: {volume} against {}",
        8000.0 + filled
    );
}

/// The faces of `shape` whose surface `keep` accepts.
fn faces_on(
    model: &Model,
    shape: &Shape,
    keep: impl Fn(&ogeom::geom::SurfaceGeometry) -> bool,
) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            keep(model.geometry().surface(data.surface).unwrap())
        })
        .collect()
}

/// The planar face of `shape` at height `z` facing up.
fn lid_at(model: &Model, shape: &Shape, z: f64) -> Shape {
    let found: Vec<Shape> = faces_on(model, shape, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::Plane(_))
    })
    .into_iter()
    .filter(|f| {
        let (p, n) = ogeom::algo::face_normal(model, f, T).unwrap();
        (p.z - z).abs() < 1e-9 && n.z > 0.5
    })
    .collect();
    assert_eq!(found.len(), 1, "one lid at z = {z}");
    found[0].clone()
}

/// Whether two faces share an edge.
fn share_an_edge(model: &Model, a: &Shape, b: &Shape) -> bool {
    let theirs = explore_unique(model, b, ShapeType::Edge).unwrap();
    explore_unique(model, a, ShapeType::Edge)
        .unwrap()
        .iter()
        .any(|e| theirs.iter().any(|t| t.node() == e.node()))
}

/// The blend's joins with other faces of `shape`, as `analyse_blend`
/// measures them (a closed band also meets itself across its seam).
fn joins(model: &Model, shape: &Shape, blend: &Shape) -> Vec<ogeom::fillet::BlendContact> {
    ogeom::fillet::analyse_blend(model, shape, blend, 15, T)
        .unwrap()
        .into_iter()
        .filter(|c| c.neighbour.node() != blend.node())
        .collect()
}

/// Every join of a closed-form round: on both faces, tangent within
/// 1e-5 rad.
fn assert_exactly_tangent(model: &Model, shape: &Shape, round: &Shape) {
    let found = joins(model, shape, round);
    assert!(found.len() >= 2, "the round meets both faces: {found:?}");
    for c in &found {
        assert!(c.gap < 1e-9, "a join stands off its faces: {c:?}");
        assert!(c.tangency_error < 1e-5, "a join is not tangent: {c:?}");
    }
}

/// The circular edge of `shape` at height `z`.
fn rim_at(model: &Model, shape: &Shape, z: f64) -> Shape {
    let found: Vec<Shape> = explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let data = model.node(e).unwrap().data().as_edge().unwrap();
            let Some(ogeom::topo::EdgeRepr::Curve3d { curve, .. }) = data.curve3d() else {
                return false;
            };
            matches!(
                model.geometry().curve(*curve),
                Some(ogeom::geom::Curve::Circle(c)) if (c.circle().centre().z - z).abs() < 1e-9
            )
        })
        .collect();
    assert_eq!(found.len(), 1, "one rim at z = {z}");
    found[0].clone()
}

fn exact_volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, ogeom::mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// A square corner of side `r` less the quarter disc a ball of radius `r`
/// leaves in it, turned about an axis `reach` from the corner, the corner
/// pointing away from the axis (`sign` 1) or toward it (`sign` -1): the
/// ring a rim blend takes off or fills, by Pappus.
fn rim_ring(reach: f64, r: f64, sign: f64) -> f64 {
    let pi = core::f64::consts::PI;
    let area = r * r * (1.0 - pi / 4.0);
    // The centroid's distance from the square's corner, along each side.
    let inset = r * (10.0 - 3.0 * pi) / (3.0 * (4.0 - pi));
    2.0 * pi * sign.mul_add(-inset, reach) * area
}

/// A cylinder boss fused on a block: its wall meets the block's top along
/// a circle, a concave corner. Blending the wall and the top rolls the ball
/// round the boss's foot and fills the corner with a torus.
#[test]
fn a_face_blend_fills_a_boss_s_foot_with_a_torus() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 40.0, 10.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(20.0, 20.0, 5.0), Direction::Z, Direction::X, T).unwrap();
    let post = ogeom::algo::make_cylinder(&mut model, at, 8.0, 15.0, T)
        .unwrap()
        .shape;
    let bossed = ogeom::boolean::fuse(&mut model, &block, &post, T)
        .unwrap()
        .shape;
    let walls = faces_on(&model, &bossed, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::Cylinder(_))
    });
    assert_eq!(walls.len(), 1);
    let top = lid_at(&model, &bossed, 10.0);

    let r = 2.0;
    let blended = ogeom::fillet::blend_faces(&mut model, &bossed, &walls[0], &top, r, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &blended, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let tori: Vec<(Shape, ogeom::math::Torus)> = faces_on(&model, &blended, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::Torus(_))
    })
    .into_iter()
    .map(|f| {
        let data = model.node(&f).unwrap().data().as_face().unwrap();
        let Some(ogeom::geom::SurfaceGeometry::Torus(t)) = model.geometry().surface(data.surface)
        else {
            unreachable!()
        };
        let torus = t.torus();
        (f, torus)
    })
    .collect();
    assert_eq!(tori.len(), 1, "the round is one torus");
    let (round, torus) = &tori[0];
    assert!((torus.minor_radius() - r).abs() < 1e-12);
    assert!(
        (torus.major_radius() - 10.0).abs() < 1e-12,
        "the ball's centre runs at 8 + 2"
    );
    assert_exactly_tangent(&model, &blended, round);

    // The block, the boss above it, and the ring the ball fills, its
    // corner pointing away from the axis at the wall's radius.
    let pi = core::f64::consts::PI;
    let want = 16000.0 + pi * 64.0 * 10.0 + rim_ring(8.0, r, -1.0);
    let got = exact_volume(&model, &blended);
    assert!(
        (got - want).abs() < want * 1e-9,
        "the boss's foot is filled: {got} against {want}"
    );
}

/// A drum whose top rim is chamfered: its wall and its top share no edge,
/// the chamfer's cone standing between them. Blending the wall and the top
/// rolls the ball inside the drum's convex corner, and a radius that clears
/// the chamfer leaves the drum a plain rim round: the chamfer is gone and
/// the volume is the sharp drum's less the rim ring.
#[test]
fn a_face_blend_rounds_across_a_chamfer_between_the_faces() {
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 10.0, 20.0, T)
        .unwrap()
        .shape;
    let rim = rim_at(&model, &drum, 20.0);
    let chamfered = ogeom::fillet::chamfer_edge(&mut model, &drum, &rim, 1.0, T)
        .unwrap()
        .shape;
    let wall = cylinder_face_of(&model, &chamfered);
    let top = lid_at(&model, &chamfered, 20.0);
    assert!(
        !share_an_edge(&model, &wall, &top),
        "the chamfer parts them"
    );

    let r = 4.0;
    let blended = ogeom::fillet::blend_faces(&mut model, &chamfered, &wall, &top, r, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &blended, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(
        faces_on(&model, &blended, |s| matches!(
            s,
            ogeom::geom::SurfaceGeometry::Cone(_)
        ))
        .is_empty(),
        "the round takes the whole chamfer"
    );
    let tori = faces_on(&model, &blended, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::Torus(_))
    });
    assert_eq!(tori.len(), 1, "the round is one torus");
    assert_exactly_tangent(&model, &blended, &tori[0]);
    let pi = core::f64::consts::PI;
    let want = pi * 100.0 * 20.0 - rim_ring(10.0, r, 1.0);
    let got = exact_volume(&model, &blended);
    assert!(
        (got - want).abs() < want * 1e-9,
        "a plain rim round: {got} against {want}"
    );
}

/// Two parallel drums fused side by side meet along two straight creases,
/// both concave. Blending one drum's wall against the other's fills both
/// creases with exact cylinders, and the fill is measured against the
/// section's closed form.
#[test]
fn a_face_blend_fills_between_two_drums_of_one_solid() {
    let (big, apart, r, length) = (5.0_f64, 6.0_f64, 1.0_f64, 20.0_f64);
    let mut model = Model::new();
    let at = |x: f64| Frame::new(Point::new(x, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let left = ogeom::algo::make_cylinder(&mut model, at(0.0), big, length, T)
        .unwrap()
        .shape;
    let right = ogeom::algo::make_cylinder(&mut model, at(apart), big, length, T)
        .unwrap()
        .shape;
    let pair = ogeom::boolean::fuse(&mut model, &left, &right, T)
        .unwrap()
        .shape;
    let walls = faces_on(&model, &pair, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::Cylinder(_))
    });
    assert_eq!(walls.len(), 2);

    let blended = ogeom::fillet::blend_faces(&mut model, &pair, &walls[0], &walls[1], r, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &blended, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let rounds: Vec<Shape> = faces_on(
        &model,
        &blended,
        |s| matches!(s, ogeom::geom::SurfaceGeometry::Cylinder(c) if (c.cylinder().radius() - r).abs() < 1e-12),
    );
    assert_eq!(rounds.len(), 2, "one exact round in each crease");
    for round in &rounds {
        // Its rails ride the drums; its ends meet the end faces square.
        let rails: Vec<_> = joins(&model, &blended, round)
            .into_iter()
            .filter(|c| c.tangency_error < 1.0)
            .collect();
        assert_eq!(rails.len(), 2, "two rails: {rails:?}");
        for c in &rails {
            assert!(c.gap < 1e-9 && c.tangency_error < 1e-5, "{c:?}");
        }
    }

    // One crease's section, halved by the plane between the axes, with the
    // right axis at (3, 0), the ball's centre at (0, h) (5 + r from both
    // axes) and the crease at (0, y). The triangle (0, 0), (3, 0), (0, h)
    // holds the half fill, the right drum's part of it (the triangle
    // (0, 0), (3, 0), (0, y) and the sector from the crease to the ball's
    // touch) and the ball's sector.
    let pi = core::f64::consts::PI;
    let half = apart / 2.0;
    let h = (big + r).mul_add(big + r, -(half * half)).sqrt();
    let y = big.mul_add(big, -(half * half)).sqrt();
    let sector = big * big * (h.atan2(half) - y.atan2(half)) / 2.0;
    let ball = r * r * half.atan2(h) / 2.0;
    let fill = 2.0 * 2.0 * (half * h / 2.0 - half * y / 2.0 - sector - ball) * length;
    let lens =
        2.0 * big * big * (half / big).acos() - half * (4.0 * big * big - apart * apart).sqrt();
    let want = (2.0 * pi * big * big - lens) * length + fill;
    let got = exact_volume(&model, &blended);
    assert!(
        (got - want).abs() < fill * 1e-6,
        "both creases filled: {got} against {want}"
    );
}

/// A boss whose wall is a B-spline surface (the same drum, restated) on a
/// block: the seat has no closed form, so the ball is marched round the
/// foot and the round is a fitted band. It must ride the wall and the top
/// within a tenth of a degree, and fill what the exact torus fills.
#[test]
fn a_face_blend_marches_between_a_plane_and_a_spline_face() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 40.0, 10.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(20.0, 20.0, 5.0), Direction::Z, Direction::X, T).unwrap();
    let post = ogeom::algo::make_cylinder(&mut model, at, 8.0, 15.0, T)
        .unwrap()
        .shape;
    let post = ogeom::algo::to_nurbs(&mut model, &post, T).unwrap().shape;
    let bossed = ogeom::boolean::fuse(&mut model, &block, &post, T)
        .unwrap()
        .shape;
    let wall: Vec<Shape> = faces_on(&model, &bossed, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::BSpline(_))
    })
    .into_iter()
    .filter(|f| ogeom::algo::face_normal(&model, f, T).unwrap().1.z.abs() < 0.5)
    .collect();
    assert_eq!(wall.len(), 1, "one spline wall");
    let top = lid_at(&model, &bossed, 10.0);
    let before = exact_volume(&model, &bossed);

    let r = 2.0;
    let blended = ogeom::fillet::blend_faces(&mut model, &bossed, &wall[0], &top, r, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &blended, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let bands: Vec<Shape> = faces_on(&model, &blended, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::BSpline(_))
    })
    .into_iter()
    .filter(|f| {
        let n = ogeom::algo::face_normal(&model, f, T).unwrap().1;
        n.z.abs() > 0.1 && n.z.abs() < 0.9
    })
    .collect();
    assert_eq!(bands.len(), 1, "one fitted band");
    // Tangent along both rails within a tenth of a degree, as analyse_blend
    // measures it, and each rail on its host: the exact drum of radius 8
    // and the plane z = 10 the spline wall and the lid restate.
    let found = joins(&model, &blended, &bands[0]);
    assert!(found.len() >= 2, "the band meets both faces: {found:?}");
    for c in &found {
        assert!(
            c.tangency_error < 0.1_f64.to_radians(),
            "a rail is not tangent: {c:?}"
        );
        let data = model.node(&c.edge).unwrap().data().as_edge().unwrap();
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            panic!("a rail has a curve");
        };
        let rail = model.geometry().curve(*curve).unwrap();
        for k in 0..=64 {
            use ogeom::geom::Curve3d as _;
            let t = (range.1 - range.0).mul_add(f64::from(k) / 64.0, range.0);
            let p = rail.point_at(t, T).unwrap();
            let off_drum = ((p.x - 20.0).hypot(p.y - 20.0) - 8.0).abs();
            let off_lid = (p.z - 10.0).abs();
            assert!(
                off_drum.min(off_lid) < 1e-4,
                "a rail leaves its host at {p:?}"
            );
        }
    }
    let filled = exact_volume(&model, &blended) - before;
    let want = rim_ring(8.0, r, -1.0);
    assert!(
        (filled - want).abs() < want * 1e-4,
        "the march fills what the torus does: {filled} against {want}"
    );
}

/// A block 40 by 40 by 10 with a drum of radius 6 fused on its top, the
/// drum's axis leaning `slant` degrees from upright toward +y: its wall
/// meets the top along an ellipse, a concave corner.
fn slanted_drum_on_block(model: &mut Model, slant: f64) -> Shape {
    let block = ogeom::algo::make_box(model, Frame::WORLD, (40.0, 40.0, 10.0), T)
        .unwrap()
        .shape;
    let (s, c) = slant.to_radians().sin_cos();
    let axis = Direction::new(Vector::new(0.0, s, c), T).unwrap();
    let at = Frame::new(Point::new(20.0, 20.0, 5.0), axis, Direction::X, T).unwrap();
    let drum = ogeom::algo::make_cylinder(model, at, 6.0, 15.0, T)
        .unwrap()
        .shape;
    ogeom::boolean::fuse(model, &block, &drum, T).unwrap().shape
}

/// The drum's wall of `shape`.
fn drum_wall(model: &Model, shape: &Shape) -> Shape {
    let walls = faces_on(model, shape, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::Cylinder(_))
    });
    assert_eq!(walls.len(), 1, "one wall");
    walls[0].clone()
}

/// The edges faces `a` and `b` share.
fn shared_edges(model: &Model, a: &Shape, b: &Shape) -> Vec<Shape> {
    let theirs = explore_unique(model, b, ShapeType::Edge).unwrap();
    explore_unique(model, a, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| theirs.iter().any(|t| t.node() == e.node()))
        .collect()
}

/// The slanted drum's foot rounded at radius 0.5, which parts its wall
/// from the block's top, then the wall and the top blended at `r`: the
/// blended solid and the volume the sharp solid had.
fn foot_blended(model: &mut Model, slant: f64, r: f64) -> (Shape, f64) {
    let sharp = slanted_drum_on_block(model, slant);
    let before = exact_volume(model, &sharp);
    let foot = shared_edges(
        model,
        &drum_wall(model, &sharp),
        &lid_at(model, &sharp, 10.0),
    );
    assert!(!foot.is_empty(), "the sharp drum stands on the top");
    let parted = ogeom::fillet::fillet_edges(model, &sharp, &foot, 0.5, T)
        .unwrap()
        .shape;
    let (wall, top) = (drum_wall(model, &parted), lid_at(model, &parted, 10.0));
    assert!(
        !share_an_edge(model, &wall, &top),
        "the small round parts them"
    );
    let blended = ogeom::fillet::blend_faces(model, &parted, &wall, &top, r, T)
        .unwrap()
        .shape;
    (blended, before)
}

/// A drum leaning 15 degrees on a block, its foot rounded small so its
/// wall and the block's top share no edge: their surfaces meet along an
/// ellipse and share no direction or axis. Blending them marches the ball
/// round the ellipse and fills the corner with a fitted band, which takes
/// the small round whole. The band rides both faces within a tenth of a
/// degree, and the fill is what rounding the sharp foot as an edge fills.
#[test]
fn a_face_blend_marches_round_a_slanted_drum_s_foot() {
    let r = 2.0;
    let mut model = Model::new();
    let (blended, before) = foot_blended(&mut model, 15.0, r);
    let diagnosis = ogeom::algo::check(&model, &blended, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let bands = faces_on(&model, &blended, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::BSpline(_))
    });
    assert_eq!(bands.len(), 1, "one fitted band, the small round taken");
    let found = joins(&model, &blended, &bands[0]);
    assert!(found.len() >= 2, "the band meets both faces: {found:?}");
    for c in &found {
        assert!(
            c.tangency_error < 0.1_f64.to_radians(),
            "a rail is not tangent: {c:?}"
        );
    }
    let filled = exact_volume(&model, &blended) - before;

    // The same corner rounded as the edge it is on the sharp solid.
    let mut other = Model::new();
    let sharp = slanted_drum_on_block(&mut other, 15.0);
    let foot = shared_edges(
        &other,
        &drum_wall(&other, &sharp),
        &lid_at(&other, &sharp, 10.0),
    );
    let rounded = ogeom::fillet::fillet_edges(&mut other, &sharp, &foot, r, T)
        .unwrap()
        .shape;
    let want = exact_volume(&other, &rounded) - exact_volume(&other, &sharp);
    assert!(
        (filled - want).abs() < want * 1e-4,
        "the blend fills what the edge round does: {filled} against {want}"
    );
}

/// The slanted drum's fill against the closed form it tends to: upright,
/// the round is a torus and the fill a ring by Pappus. The fill grows
/// with the square of the slant, so the fills at 1 and 2 degrees
/// extrapolate to the ring's: four times the first less the second, over
/// three.
#[test]
fn a_slanted_drum_s_fill_tends_to_the_torus_ring() {
    let r = 2.0;
    let ring = rim_ring(6.0, r, -1.0);
    let mut fills = Vec::new();
    for slant in [1.0, 2.0] {
        let mut model = Model::new();
        let (blended, before) = foot_blended(&mut model, slant, r);
        fills.push(exact_volume(&model, &blended) - before);
    }
    let (one, two) = (fills[0] - ring, fills[1] - ring);
    assert!(one > 0.0, "a slant widens the fill: {one}");
    assert!(
        (3.5..4.5).contains(&(two / one)),
        "the fill grows with the slant's square: {one} then {two}"
    );
    let extrapolated = 4.0f64.mul_add(fills[0], -fills[1]) / 3.0;
    assert!(
        (extrapolated - ring).abs() < ring * 1e-5,
        "the fills tend to the ring: {extrapolated} against {ring}"
    );
}

/// What the curved face blend does not build, refused by name.
#[test]
fn curved_face_blends_refuse_by_name() {
    let refusal = |r: ogeom::core::OgeomResult<ogeom::algo::Built>| match r {
        Ok(_) => panic!("refused"),
        Err(e) => e.to_string(),
    };
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 10.0, 20.0, T)
        .unwrap()
        .shape;
    let rim = rim_at(&model, &drum, 20.0);
    let chamfered = ogeom::fillet::chamfer_edge(&mut model, &drum, &rim, 1.0, T)
        .unwrap()
        .shape;
    let wall = cylinder_face_of(&model, &chamfered);
    let top = lid_at(&model, &chamfered, 20.0);
    // A ball whose round would cross the chamfer: the middle of the round
    // stands in the open the chamfer already cut.
    let said = refusal(ogeom::fillet::blend_faces(
        &mut model, &chamfered, &wall, &top, 1.5, T,
    ));
    assert!(
        said.contains("the solid beside the round is open"),
        "{said}"
    );
    // A ball smaller than the chamfer touches the two surfaces where the
    // chamfer has cut both faces away.
    let said = refusal(ogeom::fillet::blend_faces(
        &mut model, &chamfered, &wall, &top, 0.5, T,
    ));
    assert!(said.contains("does not touch both faces"), "{said}");
    // A ball wider than the drum's radius less its own reaches past the
    // axis.
    let said = refusal(ogeom::fillet::blend_faces(
        &mut model, &chamfered, &wall, &top, 6.0, T,
    ));
    assert!(said.contains("crossing its own axis"), "{said}");
    let said = refusal(ogeom::fillet::blend_faces(
        &mut model, &chamfered, &wall, &wall, 1.0, T,
    ));
    assert!(said.contains("against itself"), "{said}");
    let other = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 3.0, 5.0, T)
        .unwrap()
        .shape;
    let stranger = cylinder_face_of(&model, &other);
    let said = refusal(ogeom::fillet::blend_faces(
        &mut model, &chamfered, &stranger, &top, 1.0, T,
    ));
    assert!(said.contains("not a face of the solid"), "{said}");

    // A spline cap and a side of the block sharing no edge: no closed
    // form, and the cap's surface, continued past the face, never reaches
    // the side's plane, so there is no corner to march along.
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 40.0, 10.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(20.0, 20.0, 5.0), Direction::Z, Direction::X, T).unwrap();
    let post = ogeom::algo::make_cylinder(&mut model, at, 8.0, 15.0, T)
        .unwrap()
        .shape;
    let post = ogeom::algo::to_nurbs(&mut model, &post, T).unwrap().shape;
    let bossed = ogeom::boolean::fuse(&mut model, &block, &post, T)
        .unwrap()
        .shape;
    let cap: Vec<Shape> = faces_on(&model, &bossed, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::BSpline(_))
    })
    .into_iter()
    .filter(|f| ogeom::algo::face_normal(&model, f, T).unwrap().1.z > 0.5)
    .collect();
    assert_eq!(cap.len(), 1, "one spline cap");
    let side = faces_on(
        &model,
        &bossed,
        |s| matches!(s, ogeom::geom::SurfaceGeometry::Plane(p) if p.plane().normal().vector().x < -0.5),
    );
    let said = refusal(ogeom::fillet::blend_faces(
        &mut model, &bossed, &cap[0], &side[0], 1.0, T,
    ));
    assert!(said.contains("do not meet"), "{said}");

    // A drum leaning 15 degrees toward +y on a block whose top is split
    // along y = 29.5, past the drum's foot: the part of the top beyond the
    // line shares no edge with the wall, and the ball rolling round the
    // foot touches it only where its line of contact bulges past the line
    // on the side the drum leans over.
    let mut model = Model::new();
    let solid = slanted_drum_on_block(&mut model, 15.0);
    let top = lid_at(&model, &solid, 10.0);
    let line = ogeom::geom::Curve::Line(
        ogeom::geom::LineCurve::segment(
            Point::new(-1.0, 29.5, 10.0),
            Point::new(41.0, 29.5, 10.0),
            T,
        )
        .unwrap(),
    );
    let range = {
        use ogeom::geom::Curve3d as _;
        line.domain()
    };
    let cut = ogeom::algo::make_edge(&mut model, line, range, T)
        .unwrap()
        .shape;
    let split = ogeom::heal::split_face(
        &mut model,
        &solid,
        &top,
        &[cut],
        ogeom::heal::Projection::OnFace,
        T,
    )
    .unwrap()
    .shape;
    let beyond: Vec<Shape> = faces_on(&model, &split, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::Plane(_))
    })
    .into_iter()
    .filter(|f| {
        let (p, n) = ogeom::algo::face_normal(&model, f, T).unwrap();
        let low = ogeom::algo::shape_bounds(&model, f, T)
            .unwrap()
            .corners()
            .iter()
            .map(|q| q.y)
            .fold(f64::INFINITY, f64::min);
        (p.z - 10.0).abs() < 1e-9 && n.z > 0.5 && low > 29.0
    })
    .collect();
    assert_eq!(beyond.len(), 1, "one piece of the top past the line");
    let wall = drum_wall(&model, &split);
    assert!(!share_an_edge(&model, &wall, &beyond[0]));
    let said = refusal(ogeom::fillet::blend_faces(
        &mut model, &split, &wall, &beyond[0], 2.0, T,
    ));
    assert!(said.contains("part of its seat only"), "{said}");
}

/// The corner where three blends meet. Three edges of a box are filleted
/// in sequence at one vertex, and the leftover spike is rounded by the
/// ball-and-block tool: the corner block less the ball. The result is
/// measured against a closed form derived independently, by inclusion and
/// exclusion over the corner cube: within the cube every fillet prism's removal lies inside
/// the spike's, so the removed volume is three prism runs *outside* the cube
/// plus the spike itself:
///
///   V = 10³ − 3(1 − π/4) r² (10 − r) − r³ + πr³/6
///
/// which for r = 3 is 784 + 51.75π. The blend is tangent to everything it
/// rounds by construction (each contact a chart-degenerate curve or a
/// vertex of the tool's own patch), and this test exercises the tangential
/// set-aside, the degeneracy splits, and the tolerance-carrying welds at
/// once.
#[test]
fn b2_three_fillets_and_the_corner_tool_round_the_vertex() {
    let mut model = Model::new();
    let r = 3.0;
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let mut solid = block;
    for target in [
        Point::new(10.0, 10.0, 5.0),
        Point::new(10.0, 5.0, 10.0),
        Point::new(5.0, 10.0, 10.0),
    ] {
        let edge = edge_near(&model, &solid, target);
        solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, r, T)
            .unwrap()
            .shape;
    }

    let at = |p: Point| Frame::new(p, Direction::Z, Direction::X, T).unwrap();
    let corner = Point::new(10.0 - r, 10.0 - r, 10.0 - r);
    let cblock = ogeom::algo::make_box(&mut model, at(corner), (r, r, r), T)
        .unwrap()
        .shape;
    let ball = ogeom::algo::make_sphere(&mut model, at(corner), r, T)
        .unwrap()
        .shape;
    let tool = ogeom::boolean::cut(&mut model, &cblock, &ball, T)
        .unwrap()
        .shape;
    let rounded = ogeom::boolean::cut(&mut model, &solid, &tool, T)
        .unwrap()
        .shape;

    assert!(
        ogeom::algo::check(&model, &rounded, T).unwrap().is_valid(),
        "the rounded corner is a valid solid"
    );
    let pi = core::f64::consts::PI;
    let want =
        1000.0 - 3.0 * (1.0 - pi / 4.0) * r * r * (10.0 - r) - r * r * r + pi * r * r * r / 6.0;
    let mut previous = f64::INFINITY;
    for chord in [1e-3, 1e-4] {
        let fine = ogeom::mesh::Deflection::with_chord(chord).unwrap();
        let measured = ogeom::algo::volume_properties(&model, &rounded, fine, T)
            .unwrap()
            .mass;
        let error = (measured - want).abs() / want;
        assert!(
            error < previous || error < 1e-12,
            "refining the mesh brings the measurement closer: {measured} vs {want}"
        );
        // The curved area is three band runs and the octant. The inscribed
        // deficit at chord δ runs to a few δ/r of the curved volume share.
        assert!(
            error < chord * 2.0,
            "the vertex blend against its closed form at chord {chord}: \
             {measured} vs {want}"
        );
        previous = error;
    }
}

/// `round_vertex` reproduces the closed form above: same three fillets,
/// same corner, same inclusion and exclusion reference, with the
/// ball-and-block construction in the fillet crate under its own refusals
/// instead of spelled out per call site.
#[test]
fn round_vertex_reproduces_the_b2_closed_form() {
    let mut model = Model::new();
    let r = 3.0;
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    // The corner is captured while it still exists: the fillets consume the
    // tip, and the tool reads the corner's planes from wherever the
    // vertex's *point* says they are.
    let vertex = vertex_near(&model, &block, Point::new(10.0, 10.0, 10.0));
    let mut solid = block;
    for target in [
        Point::new(10.0, 10.0, 5.0),
        Point::new(10.0, 5.0, 10.0),
        Point::new(5.0, 10.0, 10.0),
    ] {
        let edge = edge_near(&model, &solid, target);
        solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, r, T)
            .unwrap()
            .shape;
    }
    let rounded = ogeom::fillet::round_vertex(&mut model, &solid, &vertex, r, T)
        .unwrap()
        .shape;
    assert!(
        ogeom::algo::check(&model, &rounded, T).unwrap().is_valid(),
        "the rounded corner is a valid solid"
    );
    let expected = 784.0 + 51.75 * core::f64::consts::PI;
    for chord in [1e-3, 1e-4] {
        let fine = ogeom::mesh::Deflection::with_chord(chord).unwrap();
        let measured = ogeom::algo::volume_properties(&model, &rounded, fine, T)
            .unwrap()
            .mass;
        let error = (measured - expected).abs() / expected;
        assert!(
            error < chord * 2.0,
            "round_vertex against the closed form at chord {chord}: \
             {measured} vs {expected} ({error:.2e})"
        );
    }
}

/// The corner tool at every corner of the box, the three fillets in a
/// different order at each: the construction is the same whichever way the
/// corner faces and whichever edge goes first. The tool's block face meets
/// a band exactly along the arc that bounds it, and the paving must read
/// that section the same way at every corner: a hair outside the block
/// face at some and inside at others splits the band at some corners and
/// leaves it whole at the rest.
#[test]
fn round_vertex_rounds_the_corner_at_any_placement() {
    let r = 3.0;
    let expected = 784.0 + 51.75 * core::f64::consts::PI;
    let fine = ogeom::mesh::Deflection::with_chord(1e-3).unwrap();
    let mut failed: Vec<(usize, String)> = Vec::new();
    for (ci, corner) in [
        (0.0, 0.0, 0.0),
        (10.0, 0.0, 0.0),
        (0.0, 10.0, 0.0),
        (10.0, 10.0, 0.0),
        (0.0, 0.0, 10.0),
        (10.0, 0.0, 10.0),
        (0.0, 10.0, 10.0),
        (10.0, 10.0, 10.0),
    ]
    .into_iter()
    .enumerate()
    {
        let mut model = Model::new();
        let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let at = Point::new(corner.0, corner.1, corner.2);
        let vertex = vertex_near(&model, &block, at);
        // The three edges' midpoints, the order rotated by the corner.
        let mut targets = [
            Point::new(5.0, corner.1, corner.2),
            Point::new(corner.0, 5.0, corner.2),
            Point::new(corner.0, corner.1, 5.0),
        ];
        targets.rotate_left(ci % 3);
        let mut solid = block;
        for target in targets {
            let edge = edge_near(&model, &solid, target);
            solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, r, T)
                .unwrap_or_else(|e| panic!("fillet at corner {ci} near {target:?}: {e}"))
                .shape;
        }
        let outcome = ogeom::fillet::round_vertex(&mut model, &solid, &vertex, r, T)
            .map_err(|e| e.to_string())
            .and_then(|rounded| {
                if !ogeom::algo::check(&model, &rounded.shape, T)
                    .unwrap()
                    .is_valid()
                {
                    return Err("not a valid solid".to_string());
                }
                ogeom::algo::volume_properties(&model, &rounded.shape, fine, T)
                    .map(|p| p.mass)
                    .map_err(|e| e.to_string())
            });
        match outcome {
            Ok(measured) => assert!(
                (measured - expected).abs() / expected < 2e-3,
                "corner {ci}: {measured} against {expected}"
            ),
            Err(e) => failed.push((ci, e)),
        }
    }
    assert!(
        failed.is_empty(),
        "placements that did not round: {failed:?}"
    );
}

/// An oblique corner: a sheared block's origin vertex, its three edges
/// filleted one after another, then the corner tool. The block is the
/// hexahedron bounded by the host planes and the three planes through the
/// ball's centre square to the edges. The patch it leaves meets its three
/// bands and three walls tangentially, and the caps at the corner are
/// consumed while the caps at the edges' far ends stand.
#[test]
fn round_vertex_rounds_an_oblique_corner() {
    let mut model = Model::new();
    let r = 2.0;
    let (a, b, c) = (
        Vector::new(20.0, 0.0, 0.0),
        Vector::new(6.0, 20.0, 0.0),
        Vector::new(3.6, 6.0, 20.0),
    );
    let block = ogeom::algo::make_parallelepiped(&mut model, Point::ORIGIN, [a, b, c], T)
        .unwrap()
        .shape;
    let vertex = vertex_near(&model, &block, Point::ORIGIN);
    let mut solid = block;
    for edge_vector in [a, b, c] {
        let edge = edge_near(&model, &solid, Point::ORIGIN + edge_vector * 0.5);
        solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, r, T)
            .unwrap()
            .shape;
    }
    let fine = ogeom::mesh::Deflection::with_chord(2e-3).unwrap();
    let before = ogeom::algo::volume_properties(&model, &solid, fine, T)
        .unwrap()
        .mass;
    let rounded = ogeom::fillet::round_vertex(&mut model, &solid, &vertex, r, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &rounded, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    // Six walls, three bands, the three far caps, and the patch.
    let faces = explore_unique(&model, &rounded, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 13, "walls, bands, far caps and the patch");
    let patch = faces
        .iter()
        .find(|f| {
            let ogeom::topo::NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            matches!(
                model.geometry().surface(data.surface),
                Some(ogeom::geom::SurfaceGeometry::Sphere(_))
            )
        })
        .expect("the corner's spherical patch");
    let contacts = ogeom::fillet::analyse_blend(&model, &rounded, patch, 15, T).unwrap();
    assert!(!contacts.is_empty());
    for contact in &contacts {
        assert!(
            contact.gap < 1e-3 && contact.tangency_error < 5e-3,
            "the patch meets its neighbour tangentially: gap {} tangency {}",
            contact.gap,
            contact.tangency_error
        );
    }
    let after = ogeom::algo::volume_properties(&model, &rounded, fine, T)
        .unwrap()
        .mass;
    assert!(
        after < before && before - after < r * r * r,
        "the corner sheds its spike and no more: {before} -> {after}"
    );
}

/// The N-support setback at a square pyramid's apex: four planes through
/// the vertex, one ball touching all four, and the corner tool's block a
/// polyhedron of eight faces. The corner is cut first and the four edges
/// then take their flush fillets one after another, each band ending on
/// the ball's rim. The other order, four fillets and then the corner,
/// dies at the third fillet, whose predecessors crash into each other at
/// the apex. The corner's volume is measured against the closed form:
/// the block, N pyramids of height `r` over the host quads, less the
/// ball's sector, whose solid angle is the apex's angular defect.
#[test]
fn round_vertex_sets_back_a_four_edge_apex() {
    let mut model = Model::new();
    let r = 1.5;
    let base_corners = [
        Point::new(-10.0, -10.0, 0.0),
        Point::new(10.0, -10.0, 0.0),
        Point::new(10.0, 10.0, 0.0),
        Point::new(-10.0, 10.0, 0.0),
    ];
    let apex = Point::new(0.0, 0.0, 15.0);
    let base = ogeom::algo::make_polygon(&mut model, &base_corners, true, T)
        .unwrap()
        .shape;
    let tip = ogeom::algo::make_vertex(&mut model, apex).shape;
    let pyramid = ogeom::offset::make_loft(&mut model, &base, &tip, T)
        .unwrap()
        .shape;
    let vertex = vertex_near(&model, &pyramid, apex);
    let fine = ogeom::mesh::Deflection::with_chord(2e-3).unwrap();
    let volume = |model: &Model, shape: &Shape| {
        ogeom::algo::volume_properties(model, shape, fine, T)
            .unwrap()
            .mass
    };
    let before = volume(&model, &pyramid);

    let rounded = ogeom::fillet::round_vertex(&mut model, &pyramid, &vertex, r, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &rounded, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert_eq!(
        explore_unique(&model, &rounded, ShapeType::Face)
            .unwrap()
            .len(),
        10,
        "five walls, the patch and four flush ends"
    );

    // The closed form, from the pyramid's own geometry: inward normals of
    // the four lateral planes, the ball's centre on the axis a radius in
    // from each, the feet of the centre on the edges and its touch points
    // on the planes.
    let edges: [Vector; 4] =
        std::array::from_fn(|k| (base_corners[k] - apex).normalized(T).unwrap());
    let inward: [Vector; 4] = std::array::from_fn(|k| {
        let n = (base_corners[(k + 1) % 4] - base_corners[k])
            .cross(apex - base_corners[k])
            .normalized(T)
            .unwrap();
        if n.dot(Point::ORIGIN - base_corners[k]) > 0.0 {
            n
        } else {
            -n
        }
    });
    let centre = apex + Vector::new(0.0, 0.0, r / inward[0].z);
    for m in &inward {
        assert!(
            (m.dot(centre - apex) - r).abs() < 1e-9,
            "one ball touches all four"
        );
    }
    let foot = |k: usize| apex + edges[k] * (centre - apex).dot(edges[k]);
    let touch = |k: usize| centre - inward[k] * r;
    let mut block = 0.0;
    let mut defect = core::f64::consts::TAU;
    for k in 0..4 {
        let (a, b, c, d) = (apex, foot(k), touch(k), foot((k + 1) % 4));
        let area = 0.5 * ((b - a).cross(c - a).magnitude() + (c - a).cross(d - a).magnitude());
        block += r * area / 3.0;
        defect -= edges[k].dot(edges[(k + 1) % 4]).acos();
    }
    let expected = block - r * r * r * defect / 3.0;
    let after_corner = volume(&model, &rounded);
    let shed = before - after_corner;
    assert!(
        (shed - expected).abs() < expected * 5e-3,
        "the apex sheds its block less the ball's sector: {shed} vs {expected}"
    );

    // The four flush fillets after the corner, each on the edge's remaining
    // run, each shedding the same volume as the others.
    let mut solid = rounded;
    let mut shed_by_band = Vec::new();
    for (k, base_corner) in base_corners.iter().enumerate() {
        let edge = edge_near(&model, &solid, foot(k).midpoint(*base_corner));
        let was = volume(&model, &solid);
        solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, r, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom::algo::check(&model, &solid, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "after fillet {k}: {:?}",
            diagnosis.problems
        );
        shed_by_band.push(was - volume(&model, &solid));
    }
    for (k, shed) in shed_by_band.iter().enumerate() {
        assert!(
            (shed - shed_by_band[0]).abs() < shed_by_band[0] * 1e-3,
            "band {k} sheds what band 0 does: {shed} vs {}",
            shed_by_band[0]
        );
    }
    let faces = explore_unique(&model, &solid, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 10, "five walls, four bands and the patch");
    // Every blend face (the patch and the four bands) meets each of its
    // neighbours tangentially: the bands their walls and the patch, the
    // patch its four bands.
    let mut blends = 0;
    for face in &faces {
        let ogeom::topo::NodeData::Face(data) = model.node(face).unwrap().data() else {
            continue;
        };
        if matches!(
            model.geometry().surface(data.surface),
            Some(ogeom::geom::SurfaceGeometry::Plane(_))
        ) {
            continue;
        }
        blends += 1;
        let contacts = ogeom::fillet::analyse_blend(&model, &solid, face, 15, T).unwrap();
        let mut tangent = 0;
        for contact in &contacts {
            // A band's flush run-out through the base is a cut, not a join:
            // its edge lies in the base plane and meets it at an angle.
            let on_base = explore_unique(&model, &contact.edge, ShapeType::Vertex)
                .unwrap()
                .iter()
                .all(|v| {
                    model
                        .node(v)
                        .unwrap()
                        .data()
                        .as_vertex()
                        .unwrap()
                        .point
                        .z
                        .abs()
                        < 1e-6
                });
            if on_base {
                continue;
            }
            tangent += 1;
            assert!(
                contact.gap < 1e-3 && contact.tangency_error < 5e-3,
                "a blend meets its neighbour tangentially: gap {} tangency {}",
                contact.gap,
                contact.tangency_error
            );
        }
        assert!(
            tangent >= 3,
            "a band meets two walls and the patch; the patch four bands"
        );
    }
    assert_eq!(blends, 5, "the patch and four bands");
}

/// A pyramid over a polygon: its lateral planes' inward normals, for the
/// closed-form checks below.
fn inward_normals(base: &[Point], apex: Point) -> Vec<Vector> {
    let centre = base
        .iter()
        .fold(Vector::ZERO, |acc, p| acc + (*p - Point::ORIGIN))
        * (1.0 / f64::from(u32::try_from(base.len()).unwrap_or(u32::MAX)));
    let inside = Point::ORIGIN + centre;
    (0..base.len())
        .map(|k| {
            let n = (base[(k + 1) % base.len()] - base[k])
                .cross(apex - base[k])
                .normalized(T)
                .unwrap();
            if n.dot(inside - base[k]) > 0.0 { n } else { -n }
        })
        .collect()
}

/// The edge whose two vertices lie on the line through `a` and `b`: what
/// is left of a pyramid's edge once its apex is rounded, where a midpoint
/// search can land on a shorter edge beside the patch.
fn edge_along(model: &Model, shape: &Shape, a: Point, b: Point) -> Shape {
    let d = b - a;
    let on_line = |p: Point| {
        let w = p - a;
        (w - d * (w.dot(d) / d.dot(d))).magnitude() < 1e-6
    };
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let ends = model.children_of(e).unwrap();
            ends.len() == 2
                && ends
                    .iter()
                    .all(|v| on_line(model.node(v).unwrap().data().as_vertex().unwrap().point))
        })
        .unwrap()
}

/// A pyramid's apex rounded: the solid, the vertex, and the volume before.
fn pyramid_apex(model: &mut Model, base: &[Point], apex: Point) -> (Shape, Shape) {
    let polygon = ogeom::algo::make_polygon(model, base, true, T)
        .unwrap()
        .shape;
    let tip = ogeom::algo::make_vertex(model, apex).shape;
    let pyramid = ogeom::offset::make_loft(model, &polygon, &tip, T)
        .unwrap()
        .shape;
    let vertex = vertex_near(model, &pyramid, apex);
    (pyramid, vertex)
}

/// A sphere among a shape's faces: its centre and radius.
type Ball = (Point, f64);
/// A cylinder among a shape's faces: its axis point, direction and radius.
type Drum = (Point, Vector, f64);

/// The spheres and cylinders among a shape's faces.
fn balls_and_drums(model: &Model, shape: &Shape) -> (Vec<Ball>, Vec<Drum>) {
    let mut balls = Vec::new();
    let mut drums = Vec::new();
    for face in explore_unique(model, shape, ShapeType::Face).unwrap() {
        let ogeom::topo::NodeData::Face(data) = model.node(&face).unwrap().data() else {
            continue;
        };
        match model.geometry().surface(data.surface) {
            Some(ogeom::geom::SurfaceGeometry::Sphere(s)) => {
                balls.push((s.sphere().frame().origin(), s.sphere().radius()));
            }
            Some(ogeom::geom::SurfaceGeometry::Cylinder(c)) => {
                drums.push((
                    c.cylinder().frame().origin(),
                    c.cylinder().frame().z().vector(),
                    c.cylinder().radius(),
                ));
            }
            _ => {}
        }
    }
    (balls, drums)
}

/// A rectangular pyramid's apex: four planes, two slopes, and no ball a
/// radius in from all four at once. The region the ball's centre may
/// occupy has two tip vertices, each a radius in from three of the planes,
/// joined along the two long slopes. The rounded corner is a sphere at
/// each and a cylinder between them, the exact envelope of the rolling
/// ball. The four edges then take their flush fillets, each band ending on
/// its sphere's rim.
#[test]
fn round_vertex_rounds_an_apex_no_ball_touches() {
    let mut model = Model::new();
    let r = 1.5;
    let base = [
        Point::new(-10.0, -5.0, 0.0),
        Point::new(10.0, -5.0, 0.0),
        Point::new(10.0, 5.0, 0.0),
        Point::new(-10.0, 5.0, 0.0),
    ];
    let apex = Point::new(0.0, 0.0, 15.0);
    let (pyramid, vertex) = pyramid_apex(&mut model, &base, apex);
    let fine = ogeom::mesh::Deflection::with_chord(2e-3).unwrap();
    let volume = |model: &Model, shape: &Shape| {
        ogeom::algo::volume_properties(model, shape, fine, T)
            .unwrap()
            .mass
    };
    let before = volume(&model, &pyramid);

    let rounded = ogeom::fillet::round_vertex(&mut model, &pyramid, &vertex, r, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &rounded, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert_eq!(
        explore_unique(&model, &rounded, ShapeType::Face)
            .unwrap()
            .len(),
        12,
        "five walls, two spheres, the ridge's cylinder and four flush ends"
    );
    // Each sphere sits a radius in from three of the four planes and
    // farther from the fourth. The cylinder runs between the two centres
    // along the two long slopes, at the same radius.
    let inward = inward_normals(&base, apex);
    let (balls, drums) = balls_and_drums(&model, &rounded);
    assert_eq!(balls.len(), 2);
    assert_eq!(drums.len(), 1);
    for (centre, radius) in &balls {
        assert!((radius - r).abs() < 1e-9);
        let distances: Vec<f64> = inward.iter().map(|m| m.dot(*centre - apex)).collect();
        let touching = distances.iter().filter(|d| (*d - r).abs() < 1e-7).count();
        assert_eq!(
            touching, 3,
            "a tip vertex touches three planes: {distances:?}"
        );
        assert!(distances.iter().all(|d| *d >= r - 1e-7));
    }
    let (origin, axis, radius) = drums[0];
    assert!((radius - r).abs() < 1e-9);
    for (centre, _) in &balls {
        let off = (*centre - origin) - axis * (*centre - origin).dot(axis);
        assert!(
            off.magnitude() < 1e-7,
            "the ridge runs through both centres"
        );
    }
    let after_corner = volume(&model, &rounded);
    assert!(after_corner < before && before - after_corner < 20.0 * r * r * r);

    // The flush fillets follow, each consuming its flush end.
    let mut solid = rounded;
    let mut last = after_corner;
    for corner in &base {
        let edge = edge_near(
            &model,
            &solid,
            Point::new(corner.x, corner.y, 0.0)
                + (apex - Point::new(corner.x, corner.y, 0.0)) * 0.5,
        );
        solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, r, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom::algo::check(&model, &solid, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        let now = volume(&model, &solid);
        assert!(now < last, "each band sheds material: {last} -> {now}");
        last = now;
    }
    let (balls, drums) = balls_and_drums(&model, &solid);
    assert_eq!(
        (balls.len(), drums.len()),
        (2, 5),
        "two spheres, the ridge and four bands"
    );
    assert_eq!(
        explore_unique(&model, &solid, ShapeType::Face)
            .unwrap()
            .len(),
        12
    );
}

/// A flat rectangular pyramid's apex, and an oblique one: the same two
/// spheres and a ridge, the edges taking their flush fillets after. At the
/// oblique apex the corner's sphere clears the fourth plane by a few
/// hundredths of a millimetre, so the envelope keeps a sliver of that
/// plane beside the patch, and the flush fillet that meets the sliver runs
/// a straight end tangent to the sphere's rim there, at twice the radius
/// as at the first.
#[test]
fn round_vertex_rounds_flat_and_oblique_apexes_no_ball_touches() {
    let r = 1.5;
    let base = [
        Point::new(-10.0, -4.0, 0.0),
        Point::new(10.0, -4.0, 0.0),
        Point::new(10.0, 4.0, 0.0),
        Point::new(-10.0, 4.0, 0.0),
    ];
    for (apex, r) in [
        (Point::new(0.0, 0.0, 8.0), r),
        (Point::new(3.0, 1.0, 15.0), r),
        (Point::new(3.0, 1.0, 15.0), 2.0 * r),
    ] {
        let mut model = Model::new();
        let (pyramid, vertex) = pyramid_apex(&mut model, &base, apex);
        let rounded = ogeom::fillet::round_vertex(&mut model, &pyramid, &vertex, r, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom::algo::check(&model, &rounded, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        let inward = inward_normals(&base, apex);
        let (balls, drums) = balls_and_drums(&model, &rounded);
        assert_eq!((balls.len(), drums.len()), (2, 1));
        for (centre, _) in &balls {
            let touching = inward
                .iter()
                .filter(|m| (m.dot(*centre - apex) - r).abs() < 1e-7)
                .count();
            assert_eq!(touching, 3);
        }
        let mut solid = rounded;
        for corner in &base {
            let edge = edge_along(&model, &solid, *corner, apex);
            solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, r, T)
                .unwrap()
                .shape;
            assert!(ogeom::algo::check(&model, &solid, T).unwrap().is_valid());
        }
        let (balls, drums) = balls_and_drums(&model, &solid);
        assert_eq!((balls.len(), drums.len()), (2, 5), "{apex:?} at {r}");
    }
}

/// An irregular pentagonal pyramid's apex: five planes, three tip vertices
/// and two ridges, one of them seven microns long, a sliver of cylinder
/// the tool keeps rather than merging into a sphere that touches none of
/// its planes exactly. The five edges then take their flush fillets.
#[test]
fn round_vertex_rounds_a_five_edged_apex_with_a_sliver_ridge() {
    let mut model = Model::new();
    let base: Vec<Point> = (0..5)
        .map(|k| {
            let a = core::f64::consts::TAU * f64::from(k) / 5.0 + 0.3;
            let radius = if k % 2 == 0 { 10.0 } else { 7.0 };
            Point::new(radius * a.cos(), radius * a.sin(), 0.0)
        })
        .collect();
    let apex = Point::new(1.0, 0.5, 14.0);
    let (pyramid, vertex) = pyramid_apex(&mut model, &base, apex);
    let rounded = ogeom::fillet::round_vertex(&mut model, &pyramid, &vertex, 1.5, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &rounded, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let inward = inward_normals(&base, apex);
    let (balls, drums) = balls_and_drums(&model, &rounded);
    assert_eq!(
        (balls.len(), drums.len()),
        (3, 2),
        "three spheres and two ridges"
    );
    for (centre, _) in &balls {
        let touching = inward
            .iter()
            .filter(|m| (m.dot(*centre - apex) - 1.5).abs() < 1e-7)
            .count();
        assert_eq!(touching, 3);
        assert!(inward.iter().all(|m| m.dot(*centre - apex) >= 1.5 - 1e-7));
    }
    assert_eq!(
        explore_unique(&model, &rounded, ShapeType::Face)
            .unwrap()
            .len(),
        16,
        "six walls, three spheres, two ridges and five flush ends"
    );
    let mut solid = rounded;
    for corner in &base {
        let edge = edge_near(&model, &solid, *corner + (apex - *corner) * 0.5);
        solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, 1.5, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom::algo::check(&model, &solid, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    }
    let (balls, drums) = balls_and_drums(&model, &solid);
    assert_eq!(
        (balls.len(), drums.len()),
        (3, 7),
        "five bands join the corner"
    );
}

/// A vertex where only two surfaces meet is no corner: a drum's rim meets
/// its own seam there, and the refusal says what a corner needs.
#[test]
fn round_vertex_refuses_a_vertex_of_two_surfaces_by_name() {
    let mut model = Model::new();
    let cyl = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
        .unwrap()
        .shape;
    let v = vertex_near(&model, &cyl, Point::new(5.0, 0.0, 10.0));
    let err = ogeom::fillet::round_vertex(&mut model, &cyl, &v, 1.0, T)
        .unwrap_err()
        .to_string();
    assert!(
        err.contains("three faces meet"),
        "the refusal names what a corner needs: {err}"
    );
}

/// A cube of side 10 in common with a drum of radius 6 standing on its
/// middle: the drum shaves the four vertical corners, and each top corner
/// is where the top, a side and the drum meet: a corner with a curved face
/// and a curved edge. Returns the part, the corner and the three edges'
/// midpoints (top and side, top and drum, side and drum).
fn shaved_cube(model: &mut Model) -> (Shape, Point, [Point; 3]) {
    let block = ogeom::algo::make_box(model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(5.0, 5.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drum = ogeom::algo::make_cylinder(model, frame, 6.0, 12.0, T)
        .unwrap()
        .shape;
    let part = ogeom::boolean::common(model, &block, &drum, T)
        .unwrap()
        .shape;
    let y = 5.0 - 11.0_f64.sqrt();
    let arc_mid = {
        // Halfway round the drum's rim between the corner and the side y = 0.
        let a = (y - 5.0).atan2(-5.0);
        let b = (-5.0_f64).atan2(y - 5.0);
        let m = f64::midpoint(a, b);
        Point::new(5.0 + 6.0 * m.cos(), 5.0 + 6.0 * m.sin(), 10.0)
    };
    (
        part,
        Point::new(0.0, y, 10.0),
        [Point::new(0.0, 5.0, 10.0), arc_mid, Point::new(0.0, y, 5.0)],
    )
}

/// The spherical faces of a shape and their centres and radii.
fn spheres(model: &Model, shape: &Shape) -> Vec<(Point, f64)> {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter_map(|f| {
            let data = model.node(&f)?.data().as_face()?;
            match model.geometry().surface(data.surface)? {
                ogeom::geom::SurfaceGeometry::Sphere(s) => {
                    Some((s.sphere().frame().origin(), s.sphere().radius()))
                }
                _ => None,
            }
        })
        .collect()
}

/// The ball a radius in from the top (z = 10), the side (x = 0) and the
/// drum (radius 6 about x = y = 5): its centre is a radius from each.
fn assert_corner_ball(centre: Point, radius: f64) {
    let from_drum = 6.0 - (centre.x - 5.0).hypot(centre.y - 5.0);
    for (host, d) in [
        ("top", 10.0 - centre.z),
        ("side", centre.x),
        ("drum", from_drum),
    ] {
        assert!(
            (d - radius).abs() < 1e-9,
            "the corner ball stands {d} from the {host}, not {radius}"
        );
    }
}

#[test]
fn a_corner_with_a_curved_face_rounds_with_one_ball() {
    let mut model = Model::new();
    let (part, corner, _) = shaved_cube(&mut model);
    let v = vertex_near(&model, &part, corner);
    let rounded = ogeom::fillet::round_vertex(&mut model, &part, &v, 1.0, T).unwrap();
    let diagnosis = ogeom::algo::check(&model, &rounded.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let balls = spheres(&model, &rounded.shape);
    assert_eq!(balls.len(), 1, "one ball rounds the corner");
    assert!((balls[0].1 - 1.0).abs() < 1e-12);
    assert_corner_ball(balls[0].0, 1.0);
    // The corner's tip is gone.
    let deflection = ogeom::mesh::Deflection::default();
    assert_eq!(
        ogeom::algo::classify_in_solid(&model, &rounded.shape, corner, deflection, T).unwrap(),
        ogeom::algo::Containment::Out
    );
}

#[test]
fn a_curved_corner_closes_the_same_way_round_either_order() {
    // Two routes to one rounded corner: all three edges in one call (the
    // corner first, the bands stopping flush against its ball), and the
    // bands one at a time with the corner tool after. The drum's ruling
    // goes first on the second route. Its band and the top's two meet
    // at the corner only through the ball.
    let mut model = Model::new();
    let (part, corner, mids) = shaved_cube(&mut model);
    let edges: Vec<Shape> = mids.iter().map(|m| edge_near(&model, &part, *m)).collect();
    let together = ogeom::fillet::fillet_edges(&mut model, &part, &edges, 1.0, T).unwrap();
    let diagnosis = ogeom::algo::check(&model, &together.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let balls = spheres(&model, &together.shape);
    assert_eq!(balls.len(), 1, "one ball rounds the corner");
    assert_corner_ball(balls[0].0, 1.0);

    let v = vertex_near(&model, &part, corner);
    let mut current = part.clone();
    for i in [2, 0, 1] {
        let live = edge_near(&model, &current, mids[i]);
        current = ogeom::fillet::fillet_edge(&mut model, &current, &live, 1.0, T)
            .unwrap()
            .shape;
    }
    let apart = ogeom::fillet::round_vertex(&mut model, &current, &v, 1.0, T).unwrap();
    let diagnosis = ogeom::algo::check(&model, &apart.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // Every face of both is analytic, so the volumes integrate exactly and
    // agree to rounding.
    let volume = |shape: &Shape| {
        ogeom::algo::volume_properties(&model, shape, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass
    };
    let (a, b) = (volume(&together.shape), volume(&apart.shape));
    assert!((a - b).abs() < 1e-9 * a, "{a} one call, {b} one at a time");
}

/// The rim arc's band first, then the straight band and the ruling's in
/// either order: the straight band's section through the rim's torus runs
/// tangent to the torus's end meridian into the corner below the rim, and
/// ends at that vertex. The corner tool then closes the same solid the
/// ruling-first order does.
#[test]
fn a_curved_corner_closes_with_the_rim_s_band_first() {
    let fine = ogeom::mesh::Deflection::with_chord(1e-3).unwrap();
    let mut volumes = Vec::new();
    for order in [[2, 0, 1], [1, 0, 2], [1, 2, 0]] {
        let mut model = Model::new();
        let (part, corner, mids) = shaved_cube(&mut model);
        let v = vertex_near(&model, &part, corner);
        let mut current = part.clone();
        for i in order {
            let live = edge_near(&model, &current, mids[i]);
            current = ogeom::fillet::fillet_edge(&mut model, &current, &live, 1.0, T)
                .unwrap_or_else(|e| panic!("{order:?} at {i}: {e}"))
                .shape;
            let diagnosis = ogeom::algo::check(&model, &current, T).unwrap();
            assert!(
                diagnosis.is_valid(),
                "{order:?} at {i}: {:?}",
                diagnosis.problems
            );
        }
        let rounded = ogeom::fillet::round_vertex(&mut model, &current, &v, 1.0, T)
            .unwrap_or_else(|e| panic!("{order:?} corner: {e}"))
            .shape;
        let diagnosis = ogeom::algo::check(&model, &rounded, T).unwrap();
        assert!(diagnosis.is_valid(), "{order:?}: {:?}", diagnosis.problems);
        volumes.push(
            ogeom::algo::volume_properties(&model, &rounded, fine, T)
                .unwrap()
                .mass,
        );
    }
    for v in &volumes {
        assert!((v - volumes[0]).abs() < 1e-7 * volumes[0], "{volumes:?}");
    }
}

/// The straight band, then the rim arc's, then the ruling's: the rim's
/// section through the straight band crosses the point on top where the
/// cap, the band and the torus all touch, and is fitted there as closely
/// as elsewhere. Three of its faces are measured from a mesh, so the
/// volume agrees with the other orders to the mesh's resolution.
#[test]
fn a_curved_corner_closes_with_the_straight_band_first() {
    let fine = ogeom::mesh::Deflection::with_chord(1e-3).unwrap();
    let mut measured = Vec::new();
    for order in [[2, 0, 1], [0, 1, 2]] {
        let mut model = Model::new();
        let (part, corner, mids) = shaved_cube(&mut model);
        let v = vertex_near(&model, &part, corner);
        let mut current = part.clone();
        for i in order {
            let live = edge_near(&model, &current, mids[i]);
            current = ogeom::fillet::fillet_edge(&mut model, &current, &live, 1.0, T)
                .unwrap_or_else(|e| panic!("{order:?} at {i}: {e}"))
                .shape;
        }
        let rounded = ogeom::fillet::round_vertex(&mut model, &current, &v, 1.0, T)
            .unwrap_or_else(|e| panic!("{order:?} corner: {e}"))
            .shape;
        let diagnosis = ogeom::algo::check(&model, &rounded, T).unwrap();
        assert!(diagnosis.is_valid(), "{order:?}: {:?}", diagnosis.problems);
        measured.push(
            ogeom::algo::volume_properties(&model, &rounded, fine, T)
                .unwrap()
                .mass,
        );
    }
    assert!(
        (measured[0] - measured[1]).abs() < 1e-4 * measured[0],
        "{measured:?}"
    );
}

#[test]
fn a_wall_meeting_a_drum_along_a_ruling_blends_exactly() {
    // The side x = 0 meets the drum along a vertical ruling: every section
    // is the same, the band is a drum of the fillet's radius, and the
    // material removed is the section's area times the edge's length.
    let mut model = Model::new();
    let (part, _, mids) = shaved_cube(&mut model);
    let edge = edge_near(&model, &part, mids[2]);
    let radius = 1.0;
    let result = ogeom::fillet::fillet_edge(&mut model, &part, &edge, radius, T).unwrap();
    let diagnosis = ogeom::algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let bands: Vec<(Point, f64)> = explore_unique(&model, &result.shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter_map(|f| {
            let data = model.node(&f)?.data().as_face()?;
            match model.geometry().surface(data.surface)? {
                ogeom::geom::SurfaceGeometry::Cylinder(c) if c.cylinder().radius() < 2.0 => {
                    Some((c.cylinder().frame().origin(), c.cylinder().radius()))
                }
                _ => None,
            }
        })
        .collect();
    assert_eq!(bands.len(), 1, "the band is one exact drum");
    assert!((bands[0].1 - radius).abs() < 1e-12);

    // The section in the plane z = const: the crease at (0, y0), the ball
    // at (1, yc) touching the side at (0, yc) and the drum where the line
    // from the drum's axis through the ball's centre meets it. The area
    // between crease, touches and the ball's arc, by Green's theorem over
    // a finely sampled boundary.
    let y0 = 5.0 - 11.0_f64.sqrt();
    // The ball's centre: x = 1 and 5 from the axis (6 less the radius).
    let yc = 5.0 - (25.0_f64 - 16.0).sqrt();
    let centre: (f64, f64) = (1.0, yc);
    let on_drum = {
        let (dx, dy) = (centre.0 - 5.0, centre.1 - 5.0);
        let d = dx.hypot(dy);
        (5.0 + dx / d * 6.0, 5.0 + dy / d * 6.0)
    };
    let mut boundary: Vec<(f64, f64)> = Vec::new();
    let steps = 20_000;
    for i in 0..steps {
        let t = f64::from(i) / f64::from(steps);
        boundary.push((0.0, y0 + (yc - y0) * t));
    }
    let arc = |from: (f64, f64), to: (f64, f64), about: (f64, f64)| {
        let a = (from.1 - about.1).atan2(from.0 - about.0);
        let mut b = (to.1 - about.1).atan2(to.0 - about.0);
        while b - a > core::f64::consts::PI {
            b -= core::f64::consts::TAU;
        }
        while a - b > core::f64::consts::PI {
            b += core::f64::consts::TAU;
        }
        let r = (from.0 - about.0).hypot(from.1 - about.1);
        (0..steps)
            .map(|i| {
                let t = a + (b - a) * f64::from(i) / f64::from(steps);
                (about.0 + r * t.cos(), about.1 + r * t.sin())
            })
            .collect::<Vec<_>>()
    };
    boundary.extend(arc((0.0, yc), on_drum, centre));
    boundary.extend(arc(on_drum, (0.0, y0), (5.0, 5.0)));
    let area = boundary
        .iter()
        .zip(boundary.iter().cycle().skip(1))
        .map(|(p, q)| p.0 * q.1 - q.0 * p.1)
        .sum::<f64>()
        .abs()
        / 2.0;
    let volume = |shape: &Shape| {
        ogeom::algo::volume_properties(&model, shape, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass
    };
    let removed = volume(&part) - volume(&result.shape);
    // The sampled boundary's chords shave the arcs by under 1e-9 of area.
    assert!(
        (removed - area * 10.0).abs() < 1e-6,
        "removed {removed}, the section's {area} over the edge's 10"
    );
}
/// A box grooved by a tilted drum,
/// whose creases are ellipse arcs cut open by the box sides and split
/// again by the cylinder's own seam.
fn grooved_block(model: &mut Model) -> Shape {
    use ogeom::math::Vector;
    let block = ogeom::algo::make_box(model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let tilt = 0.35_f64;
    let axis = Direction::new(Vector::new(0.0, tilt.cos(), -tilt.sin()), T).unwrap();
    let frame = Frame::new(Point::new(10.0, -5.0, 11.5), axis, Direction::X, T).unwrap();
    let drum = ogeom::algo::make_cylinder(model, frame, 4.0, 30.0, T)
        .unwrap()
        .shape;
    ogeom::boolean::cut(model, &block, &drum, T).unwrap().shape
}

#[test]
fn an_open_seat_runs_out_through_the_wall() {
    // The bottom crease of the grooved block is an ellipse arc that meets
    // the box wall at both ends: an open seat. The band runs on past each
    // end until the ball has left the solid, and the cut trims it against
    // the wall: material comes off, the blend rides both hosts
    // tangentially, and it ends on the wall itself, not on a cap standing
    // short of it with a sliver of sharp crease behind.
    let mut model = Model::new();
    let grooved = grooved_block(&mut model);
    let before =
        ogeom::algo::volume_properties(&model, &grooved, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    let faces_before = explore_unique(&model, &grooved, ShapeType::Face)
        .unwrap()
        .len();
    let arc = edge_near(&model, &grooved, Point::new(10.0, 14.84, 0.0));
    let built = ogeom::fillet::fillet_edge(&mut model, &grooved, &arc, 1.0, T).unwrap();
    let after =
        ogeom::algo::volume_properties(&model, &built.shape, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    let removed = before - after;
    assert!(
        removed > 1.0 && removed < before * 0.05,
        "a run-out fillet removes a sliver, not a bite: {removed}"
    );
    // One new face (the band) and no caps: both ends are the wall's.
    assert_eq!(
        explore_unique(&model, &built.shape, ShapeType::Face)
            .unwrap()
            .len(),
        faces_before + 1,
        "the band is the only face the blend adds"
    );

    // The blend face is the fitted band. Its rails ride the hosts
    // tangentially and every other edge of it lies on the wall.
    use ogeom::topo::NodeData;
    let blend = explore_unique(&model, &built.shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let NodeData::Face(d) = model.node(f).unwrap().data() else {
                return false;
            };
            matches!(
                model.geometry().surface(d.surface),
                Some(ogeom::geom::SurfaceGeometry::BSpline(_))
            )
        })
        .expect("the fitted band is a face of the result");
    let contacts = ogeom::fillet::analyse_blend(&model, &built.shape, &blend, 15, T).unwrap();
    let mut smooth = 0;
    let mut on_wall = 0;
    for c in &contacts {
        assert!(c.gap < 1e-3, "a contact stands off its edge: {}", c.gap);
        if c.tangency_error < 5e-3 {
            smooth += 1;
            continue;
        }
        let NodeData::Face(d) = model.node(&c.neighbour).unwrap().data() else {
            panic!("a neighbour is a face");
        };
        let Some(ogeom::geom::SurfaceGeometry::Plane(plane)) = model.geometry().surface(d.surface)
        else {
            panic!("a non-tangent neighbour of the band is the wall, a plane");
        };
        let wall = plane.plane();
        assert!(
            wall.normal().vector().y.abs() > 0.999
                && wall.distance_to(Point::new(0.0, 20.0, 0.0)).abs() < 1e-9,
            "the band's other edges lie on the y=20 wall"
        );
        on_wall += 1;
    }
    assert_eq!(smooth, 2, "two tangent rails: {contacts:?}");
    assert!(
        on_wall >= 2,
        "the band ends on the wall at both ends: {contacts:?}"
    );
}

/// Two straight edges of a box meeting at a corner, blended together:
/// the later seat runs on through the earlier band and the cut trims the
/// two bands against each other. Each wedge removes (1 − π/4) r² per unit
/// length. The corner cell where both wedges reach is counted once, and
/// what both remove there is the cell outside both cylinders,
/// r³ (5/3 − π/2). One edge at a time stops flush instead, and keeps a
/// cap at the corner, the state the corner tool is built for.
#[test]
fn two_blends_meeting_at_a_corner_trim_each_other() {
    let (l, r) = (20.0_f64, 2.0_f64);
    let pi = core::f64::consts::PI;
    let want =
        l * l * 10.0 - (2.0 * (1.0 - pi / 4.0) * r * r * l - r * r * r * (5.0 / 3.0 - pi / 2.0));

    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (l, l, 10.0), T)
        .unwrap()
        .shape;
    let a = edge_near(&model, &block, Point::new(10.0, 20.0, 10.0));
    let b = edge_near(&model, &block, Point::new(20.0, 10.0, 10.0));
    let met = ogeom::fillet::fillet_edges(&mut model, &block, &[a.clone(), b.clone()], r, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &met, T).unwrap().is_valid());
    assert_eq!(
        explore_unique(&model, &met, ShapeType::Face).unwrap().len(),
        8,
        "six walls and two bands, no cap between them"
    );
    let mut previous = f64::INFINITY;
    for chord in [1e-3, 1e-4] {
        let fine = ogeom::mesh::Deflection::with_chord(chord).unwrap();
        let measured = ogeom::algo::volume_properties(&model, &met, fine, T)
            .unwrap()
            .mass;
        let error = (measured - want).abs() / want;
        assert!(
            error < previous || error < 1e-12,
            "refining brings it closer: {measured} vs {want}"
        );
        assert!(
            error < chord * 2.0,
            "two meeting blends against the closed form at chord {chord}: {measured} vs {want}"
        );
        previous = error;
    }

    // One edge at a time: the second stops flush at the first band, and
    // the corner cell keeps the material the meeting would have rounded.
    let first = ogeom::fillet::fillet_edge(&mut model, &block, &a, r, T).unwrap();
    let b_again = edge_near(&model, &first.shape, Point::new(20.0, 10.0, 10.0));
    let flush = ogeom::fillet::fillet_edge(&mut model, &first.shape, &b_again, r, T)
        .unwrap()
        .shape;
    assert_eq!(
        explore_unique(&model, &flush, ShapeType::Face)
            .unwrap()
            .len(),
        9,
        "flush: the second wedge's cap stands at the first band"
    );
    let fine = ogeom::mesh::Deflection::with_chord(1e-3).unwrap();
    let flush_volume = ogeom::algo::volume_properties(&model, &flush, fine, T)
        .unwrap()
        .mass;
    assert!(
        flush_volume > want + 0.5,
        "the flush corner keeps material the meeting removes"
    );
}

/// An L: a 2 mm cube with the quarter above `z = 1` and beyond `x = 1`
/// taken out, so its re-entrant edge runs along `y` at `(1, 1)`.
fn l_bracket(model: &mut Model) -> Shape {
    let block = ogeom::algo::make_box(model, Frame::WORLD, (2.0, 2.0, 2.0), T)
        .unwrap()
        .shape;
    let seat = Frame::new(Point::new(1.0, -0.5, 1.0), Direction::Z, Direction::X, T).unwrap();
    let notch = ogeom::algo::make_box(model, seat, (2.0, 3.0, 2.0), T)
        .unwrap()
        .shape;
    ogeom::boolean::cut(model, &block, &notch, T).unwrap().shape
}

/// An L-bracket's re-entrant blend against the end face's convex blends:
/// the concave band first, then the rim of the end face (the leg's top
/// edge, the concave band's own end arc, the wall's edge) blended as one
/// tangent chain. The ball rolls along the two lines and, between them,
/// on the concave cylinder: a quarter turn of a torus, whose volume
/// Pappus gives as the cross-section's area times its centroid's path.
/// Both blends measure against their closed forms at chord 1e-4.
///
/// The convex edge alone, ending at the re-entrant vertex, needs the
/// boolean to clip a contact's overlap to the edge it runs along: the
/// band's flush end lands on the wall's plane, its side on the line of the
/// end face's own edge up the wall, a length below it.
#[test]
fn a_rim_blend_rolls_over_the_bracket_s_concave_blend() {
    let r = 0.5_f64;
    let pi = core::f64::consts::PI;
    let mut model = Model::new();
    let bracket = l_bracket(&mut model);
    let fine = ogeom::mesh::Deflection::with_chord(1e-4).unwrap();
    let volume = |model: &Model, shape: &Shape| {
        ogeom::algo::volume_properties(model, shape, fine, T)
            .unwrap()
            .mass
    };
    assert!((volume(&model, &bracket) - 6.0).abs() < 1e-3);

    // The convex edge alone, ending at the re-entrant vertex.
    let top = edge_near(&model, &bracket, Point::new(1.75, 0.0, 1.0));
    let convex = ogeom::fillet::fillet_edge(&mut model, &bracket, &top, r, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &convex, T).unwrap().is_valid());
    let want = 6.0 - (1.0 - pi / 4.0) * r * r;
    assert!(
        (volume(&model, &convex) - want).abs() < want * 2e-4,
        "the convex band ends flush at the wall: {} vs {want}",
        volume(&model, &convex)
    );

    // The re-entrant blend, then the rim over it.
    let reentrant = edge_near(&model, &bracket, Point::new(1.0, 1.0, 1.0));
    let concave = ogeom::fillet::fillet_edge(&mut model, &bracket, &reentrant, r, T)
        .unwrap()
        .shape;
    let added = volume(&model, &concave) - 6.0;
    let want_added = 2.0 * (1.0 - pi / 4.0) * r * r;
    assert!(
        (added - want_added).abs() < want_added * 2e-3,
        "the concave band fills its corner: {added} vs {want_added}"
    );
    let arc_mid = 1.0 + r - r / core::f64::consts::SQRT_2;
    let rim = [
        edge_near(&model, &concave, Point::new(1.75, 0.0, 1.0)),
        edge_near(&model, &concave, Point::new(arc_mid, 0.0, arc_mid)),
        edge_near(&model, &concave, Point::new(1.0, 0.0, 1.75)),
    ];
    let rolled = ogeom::fillet::fillet_edges(&mut model, &concave, &rim, r, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &rolled, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let removed = volume(&model, &concave) - volume(&model, &rolled);
    // Two straight runs of 1 − r each, and the torus quarter: area
    // (1 − π/4) r² about the concave axis, centroid at r (11/6 − π/2) / (1 − π/4).
    let want_removed = 2.0 * (1.0 - pi / 4.0) * r * r * (1.0 - r)
        + (pi / 2.0) * (11.0 / 6.0 - pi / 2.0) * r * r * r;
    assert!(
        (removed - want_removed).abs() < want_removed * 2e-3,
        "the rim blend rolls over the concave band: {removed} vs {want_removed}"
    );
}

/// A straight seat and a marched one meeting at two corners: the box's
/// bottom edge along the wall, and the grooved block's elliptical crease
/// that runs out through that wall across it. In either order the later
/// blend's run-out walks on under the earlier band until the ball has left
/// the material, the cut trimming the two bands against each other at both
/// corners, and the two orders land on the same solid. The wall's bottom
/// edge is two pieces either side of the scoop, on one line. Only the left
/// is asked for, and only the left is blended whichever goes first.
#[test]
fn a_marched_blend_meets_a_straight_blend_at_its_corners() {
    let fine = ogeom::mesh::Deflection::with_chord(2e-3).unwrap();
    let mut volumes = Vec::new();
    for straight_first in [true, false] {
        let mut model = Model::new();
        let grooved = grooved_block(&mut model);
        let before = ogeom::algo::volume_properties(&model, &grooved, fine, T)
            .unwrap()
            .mass;
        let faces_before = explore_unique(&model, &grooved, ShapeType::Face)
            .unwrap()
            .len();
        let wall_edge = edge_near(&model, &grooved, Point::new(2.0, 20.0, 0.0));
        let crease = edge_near(&model, &grooved, Point::new(10.0, 14.84, 0.0));
        let order = if straight_first {
            [wall_edge, crease]
        } else {
            [crease, wall_edge]
        };
        let met = ogeom::fillet::fillet_edges(&mut model, &grooved, &order, 1.0, T)
            .unwrap()
            .shape;
        assert!(ogeom::algo::check(&model, &met, T).unwrap().is_valid());
        assert_eq!(
            explore_unique(&model, &met, ShapeType::Face).unwrap().len(),
            faces_before + 2,
            "two bands, no cap between them, straight first {straight_first}"
        );
        // The right-hand piece of the wall's bottom edge stays sharp: the
        // side wall at x = 20 keeps its whole rectangle.
        let side = explore_unique(&model, &met, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                let b = ogeom::algo::shape_bounds(&model, f, T).unwrap();
                b.low().is_some_and(|lo| lo.x > 19.999) && b.high().is_some_and(|hi| hi.x < 20.001)
            })
            .expect("the x = 20 wall");
        let side_area = ogeom::algo::surface_properties(&model, &side, fine, T)
            .unwrap()
            .mass;
        assert!(
            (side_area - 200.0).abs() < 1e-3,
            "the unrequested right piece is not blended: {side_area}"
        );
        let after = ogeom::algo::volume_properties(&model, &met, fine, T)
            .unwrap()
            .mass;
        assert!(after < before && after > before * 0.9);
        volumes.push(after);
    }
    // The marched band is fitted and its volume is measured from a mesh,
    // so the two orders agree to what the mesh resolves: a few parts in a
    // million of three thousand cubic millimetres.
    assert!(
        (volumes[0] - volumes[1]).abs() < 3e-2,
        "both orders round the same material: {volumes:?}"
    );
}

#[test]
fn two_seam_split_blends_meet_cap_to_cap_in_either_order() {
    // The top crease is split by the drum's seam into two arcs sharing a
    // vertex mid-scoop. Each rounds as a capped blend ending in the arc's
    // own section plane at the seam vertex. The second blend's cap meets
    // the first's in that plane, the two bands meet along the shared arc,
    // and both caps are consumed. Whichever arc goes first, the result is
    // the same closed solid with two bands and no cap.
    let mut volumes = Vec::new();
    for order in [
        [Point::new(7.3, 7.7, 10.0), Point::new(12.7, 7.7, 10.0)],
        [Point::new(12.7, 7.7, 10.0), Point::new(7.3, 7.7, 10.0)],
    ] {
        let mut model = Model::new();
        let grooved = grooved_block(&mut model);
        let faces_before = explore_unique(&model, &grooved, ShapeType::Face)
            .unwrap()
            .len();
        let before = ogeom::algo::volume_properties(
            &model,
            &grooved,
            ogeom::mesh::Deflection::with_chord(1e-3).unwrap(),
            T,
        )
        .unwrap()
        .mass;
        let first_arc = edge_near(&model, &grooved, order[0]);
        let first = ogeom::fillet::fillet_edge(&mut model, &grooved, &first_arc, 1.0, T).unwrap();
        let second_arc = edge_near(&model, &first.shape, order[1]);
        let second =
            ogeom::fillet::fillet_edge(&mut model, &first.shape, &second_arc, 1.0, T).unwrap();
        for shell in explore_unique(&model, &second.shape, ShapeType::Shell).unwrap() {
            assert!(ogeom::algo::is_shell_closed(&model, &shell).unwrap());
        }
        assert_eq!(
            explore_unique(&model, &second.shape, ShapeType::Face)
                .unwrap()
                .len(),
            faces_before + 2,
            "two bands, no caps"
        );
        let after = ogeom::algo::volume_properties(
            &model,
            &second.shape,
            ogeom::mesh::Deflection::with_chord(1e-3).unwrap(),
            T,
        )
        .unwrap()
        .mass;
        assert!(after < before && after > before * 0.9);
        volumes.push(after);
    }
    assert!(
        (volumes[0] - volumes[1]).abs() < 1e-2,
        "the order does not change the solid: {volumes:?}"
    );
}

#[test]
fn a_seam_split_crease_arc_rounds_with_run_out_caps() {
    // The top crease is split by the cylinder's own seam into two arcs
    // sharing a mid-scoop vertex. The reconstructed loop's midpoint stands
    // in cut-away territory, so the seat is probed on the crease itself.
    // The arc then marches its seat and lands as a capped blend.
    let mut model = Model::new();
    let grooved = grooved_block(&mut model);
    let arc = edge_near(&model, &grooved, Point::new(7.3, 7.7, 10.0));
    let built = ogeom::fillet::fillet_edge(&mut model, &grooved, &arc, 1.0, T).unwrap();
    // The result meshes as one closed solid.
    let volume =
        ogeom::algo::volume_properties(&model, &built.shape, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    assert!(volume > 0.0 && volume.is_finite());
}

#[test]
fn two_disjoint_run_out_blends_coexist_on_one_solid() {
    // The first capped blend lands on the bottom crease, the second on the
    // top crease's far seam-half, nowhere near the first. Sequential
    // marched blends must not disturb each other's wounds.
    let mut model = Model::new();
    let grooved = grooved_block(&mut model);
    let first_arc = edge_near(&model, &grooved, Point::new(10.0, 14.84, 0.0));
    let first = ogeom::fillet::fillet_edge(&mut model, &grooved, &first_arc, 1.0, T).unwrap();
    let second_arc = edge_near(&model, &first.shape, Point::new(7.3, 7.7, 10.0));
    let second = ogeom::fillet::fillet_edge(&mut model, &first.shape, &second_arc, 1.0, T).unwrap();
    let volume = ogeom::algo::volume_properties(
        &model,
        &second.shape,
        ogeom::mesh::Deflection::default(),
        T,
    )
    .unwrap()
    .mass;
    assert!(volume > 0.0 && volume.is_finite());
}

fn vertex_near(model: &Model, shape: &Shape, near: Point) -> Shape {
    explore_unique(model, shape, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .min_by(|a, b| {
            let at = |v: &Shape| {
                let p = model.node(v).unwrap().data().as_vertex().unwrap().point;
                v.transform(model.datums()).unwrap().apply(p)
            };
            at(a)
                .distance(near)
                .partial_cmp(&at(b).distance(near))
                .unwrap()
        })
        .unwrap()
}

/// The bracket's re-entrant edge and the leg's front edge asked together,
/// both ways round. They meet at a corner, and neither runs on through the
/// other's band: a wedge's cut would eat a fill, and a fill cannot run on
/// through a wedge's band, so whichever is asked first takes the corner
/// and the other stops flush against its rail. The two orders therefore
/// land on different solids (which is the honest picture of a corner no
/// single ball rolls around), and each is its own closed form: the fill
/// over the length the other blend leaves it, the wedge over what is left
/// of its own edge.
///
/// The run-on asks which way a neighbour's band rounds. Without that the
/// wedge would run the whole length of the leg and out the far side of the
/// wall.
#[test]
fn a_fill_and_a_wedge_asked_together_stop_at_each_other() {
    let r = 0.5_f64;
    let pi = core::f64::consts::PI;
    // What either blend moves per unit of its length, the wedge out and the
    // fill in: the square corner less the ball's quarter.
    let per_length = (1.0 - pi / 4.0) * r * r;
    let fine = ogeom::mesh::Deflection::with_chord(1e-4).unwrap();
    let mut volumes = Vec::new();
    for fill_first in [true, false] {
        let mut model = Model::new();
        let bracket = l_bracket(&mut model);
        let fill = edge_near(&model, &bracket, Point::new(1.0, 1.0, 1.0));
        let wedge = edge_near(&model, &bracket, Point::new(1.75, 0.0, 1.0));
        let asked = if fill_first {
            [fill, wedge]
        } else {
            [wedge, fill]
        };
        let built = ogeom::fillet::fillet_edges(&mut model, &bracket, &asked, r, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom::algo::check(&model, &built, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "fill first {fill_first}: {:?}",
            diagnosis.problems
        );
        assert_eq!(
            explore_unique(&model, &built, ShapeType::Face)
                .unwrap()
                .len(),
            13,
            "fill first {fill_first}: six walls, two bands, and the ends they stop on"
        );
        // Asked first, the fill runs the leg's whole depth and the wedge
        // then has the length its rail leaves. Asked second, the fill
        // starts where the wedge's own rail crosses the re-entrant edge.
        let want = if fill_first {
            6.0 + per_length * 2.0 - per_length * (1.0 - r)
        } else {
            6.0 - per_length + per_length * (2.0 - r)
        };
        let got = ogeom::algo::volume_properties(&model, &built, fine, T)
            .unwrap()
            .mass;
        assert!(
            (got - want).abs() < 1e-3,
            "fill first {fill_first}: {got} against {want}"
        );
        volumes.push(got);
    }
    assert!(
        (volumes[0] - volumes[1]).abs() > 1e-2,
        "the orders round different corners: {volumes:?}"
    );
}

/// The blend face of a result: the one face on a surface that is neither
/// of the two hosts' kinds: a fitted band, or a torus where none was.
fn blend_face_of(model: &Model, shape: &Shape) -> Shape {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let ogeom::topo::NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            matches!(
                model.geometry().surface(data.surface),
                Some(ogeom::geom::SurfaceGeometry::BSpline(_))
            )
        })
        .expect("the rolling ball left a fitted band")
}

/// A marched blend on a host with no ruling to lean on (a cone, a sphere,
/// a torus) is valid, closed, and tangent to both hosts along its rails.
fn assert_marched_blend(model: &Model, before: &Shape, after: &Shape, what: &str) {
    let diagnosis = ogeom::algo::check(model, after, T).unwrap();
    assert!(diagnosis.is_valid(), "{what}: {:?}", diagnosis.problems);
    let mesh =
        ogeom::mesh::triangulate(model, after, ogeom::mesh::Deflection::default(), T).unwrap();
    assert!(
        mesh.is_closed(),
        "{what}: the blended solid is not watertight"
    );
    let fine = ogeom::mesh::Deflection {
        chord: 1e-2,
        ..ogeom::mesh::Deflection::default()
    };
    let v0 = ogeom::algo::volume_properties(model, before, fine, T)
        .unwrap()
        .mass;
    let v1 = ogeom::algo::volume_properties(model, after, fine, T)
        .unwrap()
        .mass;
    assert!(
        v1 < v0 && v1 > v0 * 0.9,
        "{what}: a blend shaves a little: {v0} -> {v1}"
    );
    let blend = blend_face_of(model, after);
    let contacts = ogeom::fillet::analyse_blend(model, after, &blend, 9, T).unwrap();
    let smooth = contacts.iter().filter(|c| c.tangency_error < 2e-2).count();
    assert!(
        smooth >= 2,
        "{what}: the band is tangent along both rails: {contacts:?}"
    );
}

/// A ball drilled off its centre: the seat is a fitted seam between the
/// sphere and the bore, a host whose seat has no closed form.
///
/// The bore runs across the ball, lifted off the equator, so the seat
/// keeps clear of both poles. The rim is picked off the sphere's own seam
/// meridian, which is a circle edge of its own the nearest-edge search
/// would otherwise land on.
#[test]
fn a_marched_blend_takes_a_sphere_host() {
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(-20.0, 0.0, 2.0), Direction::X, Direction::Z, T).unwrap();
    let bore = ogeom::algo::make_cylinder(&mut model, frame, 2.0, 40.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &ball, &bore, T)
        .unwrap()
        .shape;
    let (y, z): (f64, f64) = (1.4, 2.0 + 1.42);
    let rim = edge_near(
        &model,
        &drilled,
        Point::new((100.0 - y * y - z * z).sqrt(), y, z),
    );
    let blended = ogeom::fillet::fillet_edge(&mut model, &drilled, &rim, 1.0, T)
        .unwrap()
        .shape;
    assert_marched_blend(&model, &drilled, &blended, "sphere host");
}

/// A bore straight down a ball's axis: the rim is a full circle, a
/// parallel of the sphere. The revolved blend refuses it, since the sphere
/// is neither the cap nor the coaxial wall that blend is built on, and the
/// march takes it, as it takes any other circle.
#[test]
fn a_circular_rim_on_a_sphere_takes_the_march() {
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(0.0, 0.0, -20.0), Direction::Z, Direction::X, T).unwrap();
    let bore = ogeom::algo::make_cylinder(&mut model, frame, 3.0, 40.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &ball, &bore, T)
        .unwrap()
        .shape;
    let (r, top): (f64, f64) = (3.0, 91.0_f64.sqrt());
    let rim = edge_near(
        &model,
        &drilled,
        Point::new(r * 0.5_f64.cos(), r * 0.5_f64.sin(), top),
    );
    let blended = ogeom::fillet::fillet_edge(&mut model, &drilled, &rim, 1.0, T)
        .unwrap()
        .shape;
    assert_marched_blend(&model, &drilled, &blended, "circular rim on a sphere");
}

/// A bore down the ball beside its pole: the rim stays clear of the pole,
/// but the ball rolling round it touches the sphere along a rail that
/// passes over the pole, so the sphere's leg is the band from the rail to
/// the pole with the rim cut from it. The same bore turned onto another
/// axis keeps the pole out of the leg, and the two blends agree.
#[test]
fn a_rim_whose_rail_rounds_the_sphere_s_pole_blends() {
    let fine = ogeom::mesh::Deflection::with_chord(1e-3).unwrap();
    let mut volumes = Vec::new();
    for (origin, axis, x) in [
        (Point::new(2.0, 0.0, -20.0), Direction::Z, Direction::X),
        (Point::new(-20.0, 0.0, 2.0), Direction::X, Direction::Z),
    ] {
        let mut model = Model::new();
        let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 10.0, T)
            .unwrap()
            .shape;
        let frame = Frame::new(origin, axis, x, T).unwrap();
        let bore = ogeom::algo::make_cylinder(&mut model, frame, 1.5, 40.0, T)
            .unwrap()
            .shape;
        let drilled = ogeom::boolean::cut(&mut model, &ball, &bore, T)
            .unwrap()
            .shape;
        // The rim on the far side along the bore's axis, beside the pole
        // for the first placement.
        let along = axis.vector();
        let side = origin + along * 20.0;
        let across = x.vector();
        let near = side + across * 1.5;
        let reach: f64 = (100.0 - (near - Point::ORIGIN).magnitude().powi(2)).sqrt();
        let rim = edge_near(&model, &drilled, near + along * reach);
        let blended = ogeom::fillet::fillet_edge(&mut model, &drilled, &rim, 1.0, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom::algo::check(&model, &blended, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        volumes.push(
            ogeom::algo::volume_properties(&model, &blended, fine, T)
                .unwrap()
                .mass,
        );
    }
    // Each blend is fitted to its own tolerance, and the two fits of one
    // band agree only that closely: a few parts in a hundred thousand of
    // the whole ball.
    assert!(
        (volumes[0] - volumes[1]).abs() < 5e-5 * volumes[0],
        "{volumes:?}"
    );
}

/// A ring drilled through its tube: the seat runs round the drill on the
/// torus, a host whose seat has no closed form.
#[test]
fn a_marched_blend_takes_a_torus_host() {
    let mut model = Model::new();
    let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, 10.0, 3.0, T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(10.0, 0.0, -10.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, frame, 1.5, 20.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &ring, &drill, T)
        .unwrap()
        .shape;
    let rim = edge_near(
        &model,
        &drilled,
        Point::new(11.5, 0.0, (9.0_f64 - 2.25).sqrt()),
    );
    let blended = ogeom::fillet::fillet_edge(&mut model, &drilled, &rim, 0.5, T)
        .unwrap()
        .shape;
    assert_marched_blend(&model, &drilled, &blended, "torus host");
}

/// A cone drilled across its axis: the seat is a fitted seam between the
/// cone and the bore, a host whose seat has no closed form.
#[test]
fn a_marched_blend_takes_a_cone_host() {
    let mut model = Model::new();
    let cone = ogeom::algo::make_cone(&mut model, Frame::WORLD, 12.0, 6.0, 20.0, T)
        .unwrap()
        .shape;
    // The bore runs across the axis at forty degrees from the cone's own
    // seam, so the seat lies clear of it.
    let turn = 0.7_f64;
    let along =
        ogeom::math::Direction::new(ogeom::math::Vector::new(turn.cos(), turn.sin(), 0.0), T)
            .unwrap();
    let frame = Frame::new(
        Point::new(-20.0 * turn.cos(), -20.0 * turn.sin(), 10.0),
        along,
        Direction::Z,
        T,
    )
    .unwrap();
    let bore = ogeom::algo::make_cylinder(&mut model, frame, 2.5, 40.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &cone, &bore, T)
        .unwrap()
        .shape;
    // The bore leaves the cone at radius(10) = 9 on its far side, above the
    // bore's own axis.
    let rim = edge_near(
        &model,
        &drilled,
        Point::new(9.0 * turn.cos(), 9.0 * turn.sin(), 12.5),
    );
    let blended = ogeom::fillet::fillet_edge(&mut model, &drilled, &rim, 0.8, T)
        .unwrap()
        .shape;
    assert_marched_blend(&model, &drilled, &blended, "cone host");
}

/// A block whose vertical edges are rounded, its bottom rim (four lines and
/// four quarter arcs, tangent in turn) chamfered or rounded: each arc's
/// blend turns about its corner's axis through its quarter alone, meeting
/// the neighbouring straight blends at their tangent sections, and the
/// volume removed is Pappus's: the blend's cross-section along the lines,
/// and turned through a quarter about each corner's axis.
#[test]
fn blends_along_a_rounded_rim_turn_each_arc_by_its_quarter() {
    use ogeom::algo::{check, make_box, tight_bounds, volume_properties};
    use ogeom::mesh::Deflection;
    let pi = core::f64::consts::PI;
    for (r, c, round) in [
        (1.0, 0.5, false),
        (2.0, 0.5, false),
        (3.0, 1.0, false),
        (1.0, 0.5, true),
        (2.0, 0.5, true),
    ] {
        let mut model = Model::new();
        let block = make_box(&mut model, Frame::WORLD, (20.0, 10.0, 5.0), T)
            .unwrap()
            .shape;
        let edges_where = |model: &Model, shape: &Shape, pick: &dyn Fn(Point, Point) -> bool| {
            explore_unique(model, shape, ShapeType::Edge)
                .unwrap()
                .into_iter()
                .filter(|e| {
                    let b = tight_bounds(model, e, T).unwrap();
                    pick(b.low().unwrap(), b.high().unwrap())
                })
                .collect::<Vec<_>>()
        };
        let upright = edges_where(&model, &block, &|lo, hi| hi.z - lo.z > 4.0);
        let rounded = ogeom::fillet::fillet_edges(&mut model, &block, &upright, r, T)
            .unwrap()
            .shape;
        let rim = edges_where(&model, &rounded, &|lo, hi| hi.z < 1e-9 && lo.z > -1e-9);
        assert_eq!(rim.len(), 8);
        let at = format!("corner {r}, blend {c}, rounded {round}");
        let blended = if round {
            ogeom::fillet::fillet_edges(&mut model, &rounded, &rim, c, T)
        } else {
            ogeom::fillet::chamfer_edges(&mut model, &rounded, &rim, c, T)
        }
        .unwrap_or_else(|e| panic!("{at}: {e}"))
        .shape;
        let diagnosis = check(&model, &blended, T).unwrap();
        assert!(diagnosis.is_valid(), "{at}: {diagnosis}");
        // The cross-section removed and its centroid's distance in from the
        // wall: a right triangle for the chamfer, the spandrel outside a
        // quarter circle for the fillet.
        let (section, inset) = if round {
            (
                c * c * (1.0 - pi / 4.0),
                c * 3.0f64.mul_add(-pi, 10.0) / 3.0f64.mul_add(-pi, 12.0),
            )
        } else {
            (c * c / 2.0, c / 3.0)
        };
        let lines = 2.0 * ((20.0 - 2.0 * r) + (10.0 - 2.0 * r));
        let corners = 4.0 * section * (r - inset) * pi / 2.0;
        let want =
            5.0f64.mul_add(-4.0 * r * r * (1.0 - pi / 4.0), 1000.0) - lines * section - corners;
        let v = volume_properties(&model, &blended, Deflection::default(), T)
            .unwrap()
            .mass;
        assert!((v - want).abs() < want * 1e-9, "{at}: {v} against {want}");
        // Top, bottom, four walls, four corners, four straight blends and
        // four turned ones.
        assert_eq!(
            explore_unique(&model, &blended, ShapeType::Face)
                .unwrap()
                .len(),
            18,
            "{at}"
        );
    }
}

/// The cylindrical face of a shape.
fn cylinder_face_of(model: &Model, shape: &Shape) -> Shape {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let ogeom::topo::NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            matches!(
                model.geometry().surface(data.surface),
                Some(ogeom::geom::SurfaceGeometry::Cylinder(_))
            )
        })
        .expect("a cylindrical face")
}

/// The face of a shape on a sphere.
fn sphere_face_of(model: &Model, shape: &Shape) -> Shape {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let ogeom::topo::NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            matches!(
                model.geometry().surface(data.surface),
                Some(ogeom::geom::SurfaceGeometry::Sphere(_))
            )
        })
        .expect("a spherical face")
}

#[test]
fn a_fillet_reports_the_curvature_step_where_it_meets_its_planes() {
    // A rolling-ball fillet of radius r meets its two planes tangentially:
    // no angle, but the plane is flat and the cylinder bends by 1/r square
    // to the line they share. At the end caps the cylinder's direction
    // square to the arc is its ruling, flat like the cap.
    for r in [2.0, 5.0] {
        let mut model = Model::new();
        let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 12.0), T)
            .unwrap()
            .shape;
        let edge = edge_near(&model, &block, Point::new(20.0, 0.0, 12.0));
        let blended = ogeom::fillet::fillet_edge(&mut model, &block, &edge, r, T)
            .unwrap()
            .shape;
        let blend = cylinder_face_of(&model, &blended);
        let contacts = ogeom::fillet::analyse_blend(&model, &blended, &blend, 9, T).unwrap();
        assert_eq!(contacts.len(), 4, "{r}");
        let mut lines = 0;
        let mut caps = 0;
        for contact in &contacts {
            if contact.tangency_error < 1e-9 {
                lines += 1;
                assert!(
                    (contact.curvature_error - 1.0 / r).abs() < 1e-9,
                    "r {r}: a plane against a cylinder steps by 1/r: {contact:?}"
                );
            } else {
                caps += 1;
                assert!(
                    contact.curvature_error < 1e-9,
                    "r {r}: flat both ways along the ruling: {contact:?}"
                );
            }
        }
        assert_eq!((lines, caps), (2, 2), "r {r}: {contacts:?}");
    }
}

/// The Bernstein coefficients of `p x² + q x³` over `[x0, x0 + 1]`.
fn cubic_bernstein(p: f64, q: f64, x0: f64) -> [f64; 4] {
    let c0 = (p + q * x0) * x0 * x0;
    let c1 = (2.0 * p).mul_add(x0, 3.0 * q * x0 * x0);
    let c2 = p + 3.0 * q * x0;
    let c3 = q;
    [
        c0,
        c0 + c1 / 3.0,
        c0 + 2.0 * c1 / 3.0 + c2 / 3.0,
        c0 + c1 + c2 + c3,
    ]
}

/// Two cubic-by-linear patches side by side, the graph of
/// `z = f(x) (1 + y / 5)` over `[-1, 1] x [0, 1]`: `f = x² + x³` on the
/// left and `f = p x² - 2 x³` on the right, meeting along the y axis.
///
/// Both pieces vanish with their slope at `x = 0`, so the join is always
/// tangent; their second derivatives there are `2` and `2p`, so it is
/// curvature-continuous exactly when `p = 1`, while the third derivatives
/// differ either way. The shared edge carries a pcurve on the right patch
/// only if `right_pcurve` asks.
///
/// Returns the shell, the left face and the right face.
fn spline_pair(model: &mut Model, p: f64, right_pcurve: bool) -> (Shape, Shape, Shape) {
    use ogeom::algo::{attach_pcurve, make_edge_between, make_face_on, make_vertex, make_wire};
    use ogeom::geom::{BSplineCurve, BSplineSurface, Curve, Line2d, LineCurve, SurfaceGeometry};
    use ogeom::math::bspline::ControlGrid;
    use ogeom::math::{KnotVector, Point2};
    use ogeom::topo::Location;

    let left = cubic_bernstein(1.0, 1.0, -1.0);
    let right = cubic_bernstein(p, -2.0, 0.0);
    let lift = |y: f64| 1.0 + y / 5.0;
    let cubic = || KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap();
    let row = |x0: f64, b: &[f64; 4], y: f64| -> Vec<Point> {
        (0..4)
            .map(|i| Point::new(x0 + f64::from(i) / 3.0, y, b[i as usize] * lift(y)))
            .collect()
    };
    let patch = |x0: f64, b: &[f64; 4]| -> SurfaceGeometry {
        let mut points = Vec::with_capacity(8);
        for i in 0..4 {
            for y in [0.0, 1.0] {
                points.push(row(x0, b, y)[i]);
            }
        }
        BSplineSurface::new(
            cubic(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap(),
            &ControlGrid::new(points, 4, 2).unwrap(),
            T,
        )
        .unwrap()
        .into()
    };

    let a = make_vertex(model, row(-1.0, &left, 0.0)[0]).shape;
    let b = make_vertex(model, Point::new(0.0, 0.0, 0.0)).shape;
    let c = make_vertex(model, row(0.0, &right, 0.0)[3]).shape;
    let d = make_vertex(model, row(-1.0, &left, 1.0)[0]).shape;
    let e = make_vertex(model, Point::new(0.0, 1.0, 0.0)).shape;
    let f = make_vertex(model, row(0.0, &right, 1.0)[3]).shape;
    let point_of = |m: &Model, v: &Shape| m.node(v).unwrap().data().as_vertex().unwrap().point;

    let along = |m: &mut Model, x0: f64, b: &[f64; 4], y: f64, from: &Shape, to: &Shape| {
        let curve: Curve = BSplineCurve::new(cubic(), row(x0, b, y), T).unwrap().into();
        make_edge_between(m, curve, (0.0, 1.0), from, to, T)
            .unwrap()
            .shape
    };
    let bottom_left = along(model, -1.0, &left, 0.0, &a, &b);
    let top_left = along(model, -1.0, &left, 1.0, &d, &e);
    let bottom_right = along(model, 0.0, &right, 0.0, &b, &c);
    let top_right = along(model, 0.0, &right, 1.0, &e, &f);
    let straight = |m: &mut Model, from: &Shape, to: &Shape| {
        let (p0, p1) = (point_of(m, from), point_of(m, to));
        let curve: Curve = LineCurve::segment(p0, p1, T).unwrap().into();
        make_edge_between(m, curve, (0.0, p0.distance(p1)), from, to, T)
            .unwrap()
            .shape
    };
    let side_left = straight(model, &a, &d);
    let shared = straight(model, &b, &e);
    let side_right = straight(model, &c, &f);

    let left_id = model.geometry_mut().add_surface(patch(-1.0, &left));
    let right_id = model.geometry_mut().add_surface(patch(0.0, &right));
    let chart = |m: &mut Model, edge: &Shape, id, from: (f64, f64), to: (f64, f64)| {
        let line =
            Line2d::segment(Point2::new(from.0, from.1), Point2::new(to.0, to.1), T).unwrap();
        attach_pcurve(m, edge, line.into(), id, Location::identity(), (0.0, 1.0)).unwrap();
    };
    chart(model, &bottom_left, left_id, (0.0, 0.0), (1.0, 0.0));
    chart(model, &top_left, left_id, (0.0, 1.0), (1.0, 1.0));
    chart(model, &side_left, left_id, (0.0, 0.0), (0.0, 1.0));
    chart(model, &shared, left_id, (1.0, 0.0), (1.0, 1.0));
    chart(model, &bottom_right, right_id, (0.0, 0.0), (1.0, 0.0));
    chart(model, &top_right, right_id, (0.0, 1.0), (1.0, 1.0));
    chart(model, &side_right, right_id, (1.0, 0.0), (1.0, 1.0));
    if right_pcurve {
        chart(model, &shared, right_id, (0.0, 0.0), (0.0, 1.0));
    }

    let left_wire = make_wire(
        model,
        &[
            bottom_left,
            shared.clone(),
            top_left.reversed(),
            side_left.reversed(),
        ],
        T,
    )
    .unwrap()
    .shape;
    let right_wire = make_wire(
        model,
        &[
            bottom_right,
            side_right,
            top_right.reversed(),
            shared.reversed(),
        ],
        T,
    )
    .unwrap()
    .shape;
    let left_face = make_face_on(model, left_id, &[left_wire], T).unwrap().shape;
    let right_face = make_face_on(model, right_id, &[right_wire], T)
        .unwrap()
        .shape;
    let shell = ogeom::algo::make_shell(model, &[left_face.clone(), right_face.clone()])
        .unwrap()
        .shape;
    (shell, left_face, right_face)
}

#[test]
fn two_spline_patches_built_curvature_continuous_report_no_curvature_step() {
    let mut model = Model::new();
    let (shell, left, _) = spline_pair(&mut model, 1.0, true);
    let contacts = ogeom::fillet::analyse_blend(&model, &shell, &left, 15, T).unwrap();
    assert_eq!(contacts.len(), 1, "{contacts:?}");
    let join = &contacts[0];
    assert!(join.gap < 1e-9, "{join:?}");
    assert!(join.tangency_error < 1e-9, "{join:?}");
    assert!(join.curvature_error < 1e-6, "{join:?}");
}

#[test]
fn two_spline_patches_tangent_but_bent_apart_report_their_step() {
    // With p = 3 the right piece bends by z_xx = 6 (1 + y / 5) at the join
    // against the left's 2 (1 + y / 5), the normal there straight up: the
    // step is 4 (1 + y / 5), largest at the station y = 1.
    let mut model = Model::new();
    let (shell, left, _) = spline_pair(&mut model, 3.0, true);
    let contacts = ogeom::fillet::analyse_blend(&model, &shell, &left, 15, T).unwrap();
    assert_eq!(contacts.len(), 1, "{contacts:?}");
    let join = &contacts[0];
    assert!(join.tangency_error < 1e-9, "{join:?}");
    assert!((join.curvature_error - 4.8).abs() < 1e-9, "{join:?}");
}

#[test]
fn a_join_with_no_chart_on_one_side_reports_no_curvature_it_cannot_measure() {
    let mut model = Model::new();
    let (shell, left, _) = spline_pair(&mut model, 1.0, false);
    let contacts = ogeom::fillet::analyse_blend(&model, &shell, &left, 15, T).unwrap();
    assert_eq!(contacts.len(), 1, "{contacts:?}");
    let join = &contacts[0];
    assert!(join.curvature_error.is_infinite(), "{join:?}");
    assert!(join.tangency_error.is_infinite(), "{join:?}");
    assert_eq!(join.stations, 0);
}

#[test]
fn spline_face_curvature_samples_match_the_graph_they_draw() {
    // The left patch is the graph z = (x² + x³)(1 + y / 5). Against its
    // upward normal a graph has Gaussian curvature
    // (z_xx z_yy - z_xy²) / w⁴ and mean curvature
    // ((1 + z_y²) z_xx - 2 z_x z_y z_xy + (1 + z_x²) z_yy) / (2 w³),
    // with w² = 1 + z_x² + z_y².
    let mut model = Model::new();
    let (_, left, _) = spline_pair(&mut model, 1.0, true);
    let samples = ogeom::fillet::face_curvature_samples(&model, &left, 6, T).unwrap();
    assert_eq!(samples.len(), 36, "every cell centre is inside a rectangle");
    for (at, c) in &samples {
        let (x, y) = (at.x, at.y);
        assert!(
            (-1.0..=0.0).contains(&x) && (0.0..=1.0).contains(&y),
            "{at:?}"
        );
        let (f, f1, f2) = (x * x + x * x * x, 2.0 * x + 3.0 * x * x, 2.0 + 6.0 * x);
        let g = 1.0 + y / 5.0;
        assert!((at.z - f * g).abs() < 1e-12, "on the graph: {at:?}");
        let (zx, zy, zxx, zxy, zyy) = (f1 * g, f / 5.0, f2 * g, f1 / 5.0, 0.0);
        let w2 = 1.0 + zx * zx + zy * zy;
        let gaussian = (zxx * zyy - zxy * zxy) / (w2 * w2);
        let mean = ((1.0 + zy * zy) * zxx - 2.0 * zx * zy * zxy + (1.0 + zx * zx) * zyy)
            / (2.0 * w2 * w2.sqrt());
        assert!((c.gaussian() - gaussian).abs() < 1e-9, "{at:?}: {c:?}");
        assert!((c.mean() - mean).abs() < 1e-9, "{at:?}: {c:?}");
        assert!(c.normal.vector().z > 0.0, "{c:?}");
    }
}

#[test]
fn a_fillet_samples_as_a_quarter_cylinder() {
    // The fillet of radius 2 on the box's edge along x at y = 0, z = 12
    // has its axis at y = 2, z = 10: every sample stands 2 off it, inside
    // the quarter the trim keeps, bending by -1/2 round (away from the
    // outward normal) and not at all along.
    let r = 2.0;
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 12.0), T)
        .unwrap()
        .shape;
    let edge = edge_near(&model, &block, Point::new(20.0, 0.0, 12.0));
    let blended = ogeom::fillet::fillet_edge(&mut model, &block, &edge, r, T)
        .unwrap()
        .shape;
    let blend = cylinder_face_of(&model, &blended);
    let samples = ogeom::fillet::face_curvature_samples(&model, &blend, 8, T).unwrap();
    assert_eq!(samples.len(), 64);
    for (at, c) in &samples {
        let radial = Vector::new(0.0, at.y - 2.0, at.z - 10.0);
        assert!((radial.magnitude() - r).abs() < 1e-9, "{at:?}");
        assert!(at.y <= 2.0 + 1e-9 && at.z >= 10.0 - 1e-9, "{at:?}");
        assert!((0.0..=40.0).contains(&at.x), "{at:?}");
        assert!(
            c.max.abs() < 1e-9 && (c.min + 1.0 / r).abs() < 1e-9,
            "{c:?}"
        );
        assert!(
            (c.normal.vector().dot(radial) / r - 1.0).abs() < 1e-9,
            "{c:?}"
        );
        assert!(c.max_direction.vector().x.abs() > 1.0 - 1e-9, "{c:?}");
    }
}

#[test]
fn a_drilled_face_gives_no_samples_in_its_hole() {
    let mut model = Model::new();
    let plate = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 5.0), T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(10.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let bore = ogeom::algo::make_cylinder(&mut model, frame, 3.0, 7.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &plate, &bore, T)
        .unwrap()
        .shape;
    let top = planar_face_at(&model, &drilled, Point::new(1.0, 1.0, 5.0));
    let samples = ogeom::fillet::face_curvature_samples(&model, &top, 10, T).unwrap();
    // A 10 by 10 grid of 2 mm cells: the centres within 3 of (10, 10) are
    // the four at (9 or 11, 9 or 11), sqrt(2) off; the next ring out, at
    // (7 or 13, 9 or 11) and the like, stands sqrt(10) off, past the hole.
    assert_eq!(samples.len(), 100 - 4);
    for (at, c) in &samples {
        assert!((at.z - 5.0).abs() < 1e-12, "{at:?}");
        assert!(at.distance(Point::new(10.0, 10.0, 5.0)) > 3.0, "{at:?}");
        assert!(c.max.abs() < 1e-12 && c.min.abs() < 1e-12, "{c:?}");
    }
}

#[test]
fn a_cavity_curves_toward_its_outward_normal_and_a_scaled_ball_bends_less() {
    let r = 4.0;
    let centre = Point::new(10.0, 10.0, 10.0);
    let mut model = Model::new();
    let cube = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let at_centre = Frame::new(centre, Direction::Z, Direction::X, T).unwrap();
    let ball = ogeom::algo::make_sphere(&mut model, at_centre, r, T)
        .unwrap()
        .shape;
    let hollow = ogeom::boolean::cut(&mut model, &cube, &ball, T)
        .unwrap()
        .shape;
    // The cavity's wall is the ball's face turned round: its outward normal
    // points into the hole, and the wall bends toward it.
    let cavity = sphere_face_of(&model, &hollow);
    let samples = ogeom::fillet::face_curvature_samples(&model, &cavity, 8, T).unwrap();
    assert!(samples.len() > 40, "{}", samples.len());
    for (at, c) in &samples {
        assert!((at.distance(centre) - r).abs() < 1e-9, "{at:?}");
        assert!((c.mean() - 1.0 / r).abs() < 1e-9, "{c:?}");
        assert!((c.gaussian() - 1.0 / (r * r)).abs() < 1e-9, "{c:?}");
        // The principal pair splits from the mean by the square root of a
        // rounding-sized difference, so it holds to that root of it.
        assert!(
            (c.max - 1.0 / r).abs() < 1e-7 && (c.min - 1.0 / r).abs() < 1e-7,
            "{c:?}"
        );
        assert!(c.normal.vector().dot(centre - *at) > 0.0, "{c:?}");
    }

    // The same ball doubled about a far point: its samples land on the
    // doubled sphere and bend by half as much, away from an outward normal.
    let doubling = ogeom::math::Transform::scaling(Point::new(-5.0, 2.0, 0.0), 2.0, T).unwrap();
    let grown = ogeom::algo::transformed(&mut model, &ball, doubling)
        .unwrap()
        .shape;
    let new_centre = doubling.apply(centre);
    let face = sphere_face_of(&model, &grown);
    let samples = ogeom::fillet::face_curvature_samples(&model, &face, 8, T).unwrap();
    assert!(samples.len() > 40, "{}", samples.len());
    for (at, c) in &samples {
        assert!((at.distance(new_centre) - 2.0 * r).abs() < 1e-9, "{at:?}");
        let k = -1.0 / (2.0 * r);
        assert!((c.mean() - k).abs() < 1e-9, "{c:?}");
        assert!((c.gaussian() - k * k).abs() < 1e-9, "{c:?}");
        assert!(
            (c.max - k).abs() < 1e-7 && (c.min - k).abs() < 1e-7,
            "{c:?}"
        );
        assert!(c.normal.vector().dot(*at - new_centre) > 0.0, "{c:?}");
    }
}

#[test]
fn face_curvature_samples_refuses_no_grid_and_a_shape_that_is_not_a_face() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (4.0, 3.0, 2.0), T)
        .unwrap()
        .shape;
    let face = planar_face_at(&model, &block, Point::new(1.0, 1.0, 2.0));
    assert!(matches!(
        ogeom::fillet::face_curvature_samples(&model, &face, 0, T),
        Err(ogeom::core::OgeomError::Construction(_))
    ));
    let edge = edge_near(&model, &block, Point::new(2.0, 0.0, 2.0));
    assert!(matches!(
        ogeom::fillet::face_curvature_samples(&model, &edge, 4, T),
        Err(ogeom::core::OgeomError::Construction(_))
    ));
}
