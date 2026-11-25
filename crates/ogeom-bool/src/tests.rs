use super::*;
use ogeom_algo::{check, make_box, make_cylinder, volume_properties};
use ogeom_math::{Direction, Frame};
use ogeom_mesh::Deflection;

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

#[test]
fn a_hug_measures_the_curve_rather_than_its_sampling_chords() {
    let circle = ogeom_math::Circle::new(Frame::WORLD, 100.0, T).unwrap();
    let curve = Curve::Circle(ogeom_geom::CircleCurve::new(circle));
    for range in [(0.0, 2.0 * PI), (-2.0 * PI, 2.0 * PI), (2.0 * PI, 4.0 * PI)] {
        for angle in [
            0.0_f64,
            PI,
            1.5 * PI,
            2.0 * PI,
            2.5 * PI,
            4.0 * PI,
            -PI,
            1.123456,
            -1e-8,
            1e-8,
            2.0 * PI - 1e-8,
            2.0 * PI + 1e-8,
        ] {
            for offset in [-0.123, 0.0, 0.123] {
                let p = Point::new(
                    (100.0 + offset) * angle.cos(),
                    (100.0 + offset) * angle.sin(),
                    0.0,
                );
                let distance = hug_distance(&curve, range, p, T).unwrap();
                assert!(
                    (distance - offset.abs()).abs() < 1e-8,
                    "{range:?} {angle} {offset}: {distance}"
                );
            }
        }
    }
    let p = curve.point_at(0.8, T).unwrap();
    let expected = p.distance(curve.point_at(0.7, T).unwrap());
    let distance = hug_distance(&curve, (0.4, 0.7), p, T).unwrap();
    assert!((distance - expected).abs() < 1e-10);
}

/// One cylinder, one cone and one torus, each spelled on frames turned
/// about the axis, slid along it (the cylinder's and the cone's), with
/// the axis reversed and on
/// either handedness: every pair's chart map carries a point of the
/// owner's chart to the target's point at the same place in space.
/// A frame of the other handedness turns its angle the other way
/// round, which the map must mirror rather than turn.
#[test]
fn a_coaxial_chart_map_lands_on_the_same_point_across_handedness() {
    let frames = |slide: f64| -> Vec<Frame> {
        let mut out = Vec::new();
        for up in [Direction::Z, -Direction::Z] {
            let x = Direction::new(ogeom_math::Vector::new(0.6, 0.8, 0.0), T).unwrap();
            let frame = Frame::new(Point::new(0.0, 0.0, slide), up, x, T).unwrap();
            out.extend([frame, frame.mirrored()]);
        }
        out
    };
    let cylinder = |f: Frame| -> SurfaceGeometry {
        ogeom_geom::CylinderSurface::new(ogeom_math::Cylinder::new(f, 2.0, T).unwrap(), (-5.0, 5.0))
            .unwrap()
            .into()
    };
    let torus = |f: Frame| -> SurfaceGeometry {
        ogeom_geom::TorusSurface::new(ogeom_math::Torus::new(f, 2.0, 0.5, T).unwrap()).into()
    };
    // One cone, its radius 2 + z / 4 at a height z on the world's axis,
    // spelled from each frame's origin with the half angle signed by
    // which way the frame's axis runs.
    let cone = |f: Frame| -> SurfaceGeometry {
        let up = f.z().vector().z;
        let height = f.origin().z;
        ogeom_geom::ConeSurface::new(
            ogeom_math::Cone::new(f, 0.25_f64.mul_add(height, 2.0), up * 0.25_f64.atan(), T)
                .unwrap(),
            (-5.0, 5.0),
        )
        .unwrap()
        .into()
    };
    let kinds: [(&dyn Fn(Frame) -> SurfaceGeometry, f64); 3] =
        [(&cylinder, 1.5), (&cone, 1.5), (&torus, 0.0)];
    for (make, slide) in kinds {
        for target_frame in [Frame::WORLD, Frame::WORLD.mirrored()] {
            let target = make(target_frame);
            for owner_frame in frames(slide) {
                let owner = make(owner_frame);
                let map = coaxial_chart_map(&owner, &target, 0.0, T).unwrap();
                for (u, v) in [(0.3, 0.2), (2.0, -1.0), (5.0, 1.1)] {
                    let at = owner.point_at(u, v, T).unwrap();
                    let uv = map.apply(Point2::new(u, v));
                    let tau = core::f64::consts::TAU;
                    let lifted = target
                        .point_at(uv.x.rem_euclid(tau), uv.y, T)
                        .or_else(|_| target.point_at(uv.x.rem_euclid(tau), uv.y.rem_euclid(tau), T))
                        .unwrap();
                    // Two closed forms of one point: rounding only.
                    assert!(
                        at.distance(lifted) < 1e-12,
                        "{:?} owner on {:?} target, ({u}, {v}): {at:?} against {lifted:?}",
                        owner_frame.handedness(),
                        target_frame.handedness()
                    );
                }
            }
        }
    }
}

/// An L-shaped prism across `y` whose step's underside lies on
/// `z = 0` over `1 <= x <= 2`: past its edge at `x = 1` the plane runs
/// into the leg below (a concave edge), past its edge at `x = 2` out of
/// the solid (a convex one). A thin piece of a face on that plane,
/// lying past each edge, reads inside past the first, outside past the
/// second, and between both it is dust whose strands disagree.
#[test]
fn a_thin_piece_reads_the_other_solid_s_wedge_at_its_edges() {
    let mut model = Model::new();
    let profile = [
        (0.0, -1.0),
        (1.0, -1.0),
        (1.0, 0.0),
        (2.0, 0.0),
        (2.0, 1.0),
        (0.0, 1.0),
    ]
    .map(|(x, z)| Point::new(x, 0.0, z));
    let wire = ogeom_algo::make_polygon(&mut model, &profile, true, T)
        .unwrap()
        .shape;
    let plane = ogeom_geom::PlaneSurface::over(
        ogeom_math::Plane::new(Frame::new(Point::ORIGIN, Direction::Y, Direction::X, T).unwrap()),
        (-5.0, 5.0),
        (-5.0, 5.0),
    )
    .unwrap();
    let section = ogeom_algo::make_face(&mut model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let step = ogeom_algo::make_prism(
        &mut model,
        &section,
        ogeom_math::Vector::new(0.0, 1.0, 0.0),
        T,
    )
    .unwrap()
    .shape;
    let floor = make_box(
        &mut model,
        Frame::new(Point::new(-1.0, -1.0, -3.0), Direction::Z, Direction::X, T).unwrap(),
        (4.0, 3.0, 3.0),
        T,
    )
    .unwrap()
    .shape;
    let theirs = gather(&model, &step, false, None, T).unwrap();
    let ours = gather(&model, &floor, false, None, T).unwrap();
    let top = ours
        .faces
        .iter()
        .find(|f| {
            matches!(f.surface, SurfaceGeometry::Plane(_))
                && f.edges
                    .iter()
                    .all(|e| e.curve.point_at(e.crange.0, T).unwrap().z.abs() < 1e-9)
        })
        .unwrap();
    // The step's edge along `y` at `x`, on `z = 0`, as a contact on
    // the floor's top, with the way the ring walks it so the piece lies
    // on the side `x` grows towards when `beyond` is positive.
    let contact = |x: f64, beyond: f64| -> (ContactRec, (f64, f64), bool) {
        let e = theirs
            .faces
            .iter()
            .flat_map(|f| f.edges.iter())
            .find(|e| {
                [e.crange.0, e.crange.1].iter().all(|&t| {
                    let p = e.curve.point_at(t, T).unwrap();
                    (p.x - x).abs() < 1e-9 && p.z.abs() < 1e-9
                })
            })
            .unwrap();
        let ends = [e.crange.0, e.crange.1]
            .map(|t| chart_point_within(top, e.curve.point_at(t, T).unwrap(), 1e-6, T).unwrap());
        let line = ogeom_geom::Line2d::segment(ends[0], ends[1], T).unwrap();
        let prange = ogeom_geom::Curve2d::domain(&line);
        let rec = ContactRec {
            curve: (*e.curve).clone(),
            crange: e.crange,
            pcurve: line.into(),
            prange,
            node: e.node,
            tolerance: e.tolerance,
            target_from_a: true,
            target_face: 0,
            bound: ogeom_math::Aabb::EMPTY,
            crossed: false,
        };
        // The piece is on the left of its ring about the chart's
        // normal, which is `+z` or `-z`.
        let t = f64::midpoint(e.crange.0, e.crange.1);
        let along = e.curve.d1_at(t, T).unwrap();
        let uv = ends[0];
        let normal = top.surface.normal_at(uv.x, uv.y, T).unwrap().vector();
        let left = normal.cross(along);
        (rec, e.crange, left.x * beyond < 0.0)
    };
    let thin = vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1e-6),
        Point2::new(0.0, 1e-6),
    ];
    let read = |strands: &[(f64, f64)]| -> Thin {
        let mut contacts = Vec::new();
        let mut ring = Vec::new();
        for (k, &(x, beyond)) in strands.iter().enumerate() {
            let (rec, range, reversed) = contact(x, beyond);
            contacts.push(rec);
            ring.push(Traversal {
                tag: Tag::Contact { contact: k, range },
                reversed,
                nodes: None,
            });
        }
        thin_by_strands(
            top,
            &[ring],
            std::slice::from_ref(&thin),
            &contacts,
            &theirs,
            T,
        )
        .unwrap()
    };
    assert!(matches!(read(&[(1.0, -1.0)]), Thin::Reads(PieceState::In)));
    assert!(matches!(read(&[(2.0, 1.0)]), Thin::Reads(PieceState::Out)));
    assert!(matches!(read(&[(1.0, -1.0), (2.0, 1.0)]), Thin::Disagree));
    // On the step's side of either edge the piece lies on its face: a
    // coincidence the partner reads, and the strand says nothing.
    assert!(matches!(read(&[(2.0, -1.0)]), Thin::Undecided));
}

/// A quarter of a cylinder wall of radius 0.002, its edges and vertices
/// widened to 0.003, its seam turned away from the quarter, one of its
/// rulings restated a period off. Where that ruling meets the arcs, its
/// ends stand a period from theirs, and that whole turn lifts within
/// the edges' radii, yet the surface moves along it: the joint is a
/// jump of branch, not one across an apex. Gathered, every pcurve of
/// the wall stands where it stands when the ruling is not restated.
#[test]
fn a_period_jump_on_a_wall_smaller_than_its_radii_is_folded_back() {
    let wall = |restate: f64| -> (Model, Vec<(EdgeKey, Point2)>) {
        let mut model = Model::new();
        let x = Direction::new(ogeom_math::Vector::new(-1.0, -1.0, 0.0), T).unwrap();
        let whole = make_cylinder(
            &mut model,
            Frame::new(Point::ORIGIN, Direction::Z, x, T).unwrap(),
            0.002,
            0.01,
            T,
        )
        .unwrap()
        .shape;
        let quadrant = make_box(
            &mut model,
            Frame::new(Point::new(0.0, 0.0, -1.0), Direction::Z, Direction::X, T).unwrap(),
            (1.0, 1.0, 2.0),
            T,
        )
        .unwrap()
        .shape;
        let part = common(&mut model, &whole, &quadrant, T).unwrap().shape;
        let face = explore_unique(&model, &part, ShapeType::Face).unwrap()[0].clone();
        let surface = model.node(&face).unwrap().data().as_face().unwrap().surface;
        let mut restated = false;
        for edge in explore_unique(&model, &face, ShapeType::Edge).unwrap() {
            model
                .widen(&edge, ogeom_core::Tolerance::new(0.003).unwrap())
                .unwrap();
            let Some(EdgeRepr::PCurve { curve, range, .. }) = model
                .node(&edge)
                .unwrap()
                .data()
                .as_edge()
                .unwrap()
                .pcurve_for(surface, edge.location())
                .cloned()
            else {
                continue;
            };
            let stated = model.geometry().pcurve(curve).unwrap().clone();
            let (a, b) = (
                stated.point_at(range.0, T).unwrap(),
                stated.point_at(range.1, T).unwrap(),
            );
            if restate == 0.0 || restated || (a.x - b.x).abs() > 1e-9 {
                continue;
            }
            let moved = stated
                .transformed(
                    &ogeom_math::Transform2::translation(ogeom_math::Vector2::new(restate, 0.0)),
                    T,
                )
                .unwrap();
            let moved = model.geometry_mut().add_pcurve(moved);
            if let Some(node) = model.node_mut(&edge)
                && let NodeData::Edge(data) = node.data_mut()
            {
                for repr in &mut data.representations {
                    if let EdgeRepr::PCurve {
                        curve, surface: s, ..
                    } = repr
                        && *s == surface
                    {
                        *curve = moved;
                    }
                }
            }
            restated = true;
        }
        assert!(restate == 0.0 || restated);
        let gathered = gather(&model, &part, false, None, T).unwrap();
        let starts = gathered
            .faces
            .iter()
            .find(|f| matches!(f.surface, SurfaceGeometry::Cylinder(_)))
            .unwrap()
            .edges
            .iter()
            .map(|e| (e.node, e.pcurve.point_at(e.prange.0, T).unwrap()))
            .collect();
        (model, starts)
    };
    let (_, want) = wall(0.0);
    for by in [core::f64::consts::TAU, -core::f64::consts::TAU] {
        let (_, got) = wall(by);
        assert_eq!(got.len(), want.len());
        for ((_, g), (_, w)) in got.iter().zip(&want) {
            assert!(
                g.distance(*w) < 1e-9,
                "moved {by}: {got:?} against {want:?}"
            );
        }
    }
}

/// An edge passing under a drill of radius 23.6 at right angles,
/// tangent to it a few microns off at the middle of the edge, is
/// crossed by the sections on the drill either side of the touch
/// 0.03 apart: the stretch between stays within five microns of the
/// drill. The two crossings are one junction at the touch, whether the
/// sections state a fit error of a micron's quarter or of a micron.
#[test]
fn a_graze_s_crossings_are_one_junction_however_tightly_fitted() {
    let radius = 23.6;
    let curve = Curve::Line(
        ogeom_geom::LineCurve::segment(Point::new(-10.0, 0.0, 0.0), Point::new(10.0, 0.0, 0.0), T)
            .unwrap(),
    );
    let crange = ogeom_geom::Curve3d::domain(&curve);
    let off_drill = |_: (bool, usize), p: Point| -> OgeomResult<f64> {
        Ok(p.x.hypot(p.z - radius) - radius + 4.6e-6)
    };
    let half = 0.0148;
    for honesty in [2.6e-4, 1e-3] {
        let paves: Vec<Pave> = [-half, half]
            .into_iter()
            .map(|x| Pave {
                t: x + 10.0,
                honesty,
                across: Some((false, 0)),
                spread: 0.0,
            })
            .collect();
        let clusters =
            cluster_paves(&curve, crange, 1e-7, &paves, Some(&off_drill), &[], T).unwrap();
        assert_eq!(clusters.len(), 1, "honesty {honesty}: {clusters:?}");
        assert!(clusters[0].at.x.abs() < half, "{clusters:?}");
    }
}

/// A box whose front top edge, y = 0 at z = 10, is drawn on the top
/// face `shift` inside where it runs, and states `loose`: the edge's
/// tolerance owed to the top face alone, as a chord beside a drum owes
/// its sag to the drum.
fn box_with_an_edge_drawn_off_its_top(model: &mut Model, shift: f64, loose: f64) -> Shape {
    let solid = make_box(model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let top = explore_unique(model, &solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let b = shape_bounds(model, f, T).unwrap();
            b.low().unwrap().z > 9.9
        })
        .unwrap();
    let top_surface = model.node(&top).unwrap().data().as_face().unwrap().surface;
    let edge = explore_unique(model, &top, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let b = shape_bounds(model, e, T).unwrap();
            b.high().unwrap().y.abs() < 1e-6
        })
        .unwrap();
    let (old, range) = match model
        .node(&edge)
        .unwrap()
        .data()
        .as_edge()
        .unwrap()
        .pcurve_on(top_surface)
    {
        Some(EdgeRepr::PCurve { curve, range, .. }) => (*curve, *range),
        _ => unreachable!("a box's edge has a plain pcurve on its top"),
    };
    let pcurve = model.geometry().pcurve(old).unwrap().clone();
    let (from, to) = (
        pcurve.point_at(range.0, T).unwrap(),
        pcurve.point_at(range.1, T).unwrap(),
    );
    // Square to the edge in the chart, toward the face's middle.
    let along = (to - from) / (to - from).magnitude();
    let mut across = ogeom_math::Vector2::new(-along.y, along.x);
    let middle = Point2::new(f64::midpoint(from.x, to.x), f64::midpoint(from.y, to.y));
    let lift = |p: Point2| {
        let SurfaceGeometry::Plane(plane) = model.geometry().surface(top_surface).unwrap() else {
            unreachable!("a box's top is a plane")
        };
        plane.point_at(p.x, p.y, T).unwrap()
    };
    if lift(middle + across).y < 0.0 {
        across = -across;
    }
    let line = ogeom_geom::Line2d::segment(from + across * shift, to + across * shift, T).unwrap();
    let domain = (0.0, from.distance(to));
    let id = model.geometry_mut().add_pcurve(PlanarCurve::Line(line));
    let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(|n| n.data_mut()) else {
        unreachable!("an edge")
    };
    for repr in &mut data.representations {
        if let EdgeRepr::PCurve {
            curve,
            surface,
            range,
            ..
        } = repr
            && *surface == top_surface
        {
            *curve = id;
            *range = domain;
        }
    }
    data.tolerance = data.tolerance.widen_to(loose);
    solid
}

/// Each face reads an edge as far as it stands off that face: the top
/// face draws its edge as its own pcurve on its plane, the front face
/// reads the same edge exactly. An edge stated loose that stands off
/// no face of it is read to its tolerance on both.
#[test]
fn a_face_reads_an_edge_as_far_as_it_stands_off_that_face() {
    let mut model = Model::new();
    let solid = box_with_an_edge_drawn_off_its_top(&mut model, 0.004, 0.005);
    let gathered = gather(&model, &solid, false, None, T).unwrap();
    let mut read = 0;
    for face in &gathered.faces {
        for e in &face.edges {
            if e.tolerance < 1e-3 {
                continue;
            }
            read += 1;
            let normal = face.surface.normal_at(0.0, 0.0, T).unwrap().vector();
            let (drawn, range, off) = e.drawn();
            let middle = drawn.point_at(f64::midpoint(range.0, range.1), T).unwrap();
            if normal.z.abs() > 0.5 {
                assert!((e.standoff - 0.004).abs() < 1e-6, "top: {}", e.standoff);
                assert!(e.image.is_some() && off == 0.0);
                assert!((middle.y - 0.004).abs() < 1e-6, "top drawn at {middle:?}");
            } else {
                assert!(e.standoff < 1e-9, "front: {}", e.standoff);
                assert!(e.image.is_none());
                assert!(middle.y.abs() < 1e-9, "front drawn at {middle:?}");
            }
        }
    }
    assert_eq!(read, 2);

    let mut model = Model::new();
    let solid = box_with_an_edge_drawn_off_its_top(&mut model, 1e-4, 0.005);
    let gathered = gather(&model, &solid, false, None, T).unwrap();
    for e in gathered.faces.iter().flat_map(|f| f.edges.iter()) {
        if e.tolerance >= 1e-3 {
            assert_eq!(e.standoff, e.tolerance);
            assert!(e.image.is_none());
        }
    }
}

/// The crossings either face of a loose edge finds, a lean apart along
/// it, are one cluster; tight, they stay two. A run of such crossings,
/// each within its spread of the last, makes no cluster longer than
/// two spreads.
#[test]
fn a_loose_edge_s_crossings_join_across_their_spread_and_no_further() {
    let curve = Curve::Line(
        ogeom_geom::LineCurve::segment(Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0), T)
            .unwrap(),
    );
    let crange = ogeom_geom::Curve3d::domain(&curve);
    let pave = |t: f64, spread: f64| Pave {
        t,
        honesty: 0.0,
        across: None,
        spread,
    };
    let two = [pave(1.0, 0.03), pave(1.05, 0.03)];
    let clusters = cluster_paves(&curve, crange, 0.0096, &two, None, &[], T).unwrap();
    assert_eq!(clusters.len(), 1, "{clusters:?}");
    assert!((clusters[0].span - 0.05).abs() < 1e-9, "{clusters:?}");
    let tight = [pave(1.0, 0.0), pave(1.05, 0.0)];
    assert_eq!(
        cluster_paves(&curve, crange, 0.0096, &tight, None, &[], T)
            .unwrap()
            .len(),
        2
    );
    let run: Vec<Pave> = (0..8)
        .map(|k| pave(0.5 + 0.05 * f64::from(k), 0.03))
        .collect();
    let clusters = cluster_paves(&curve, crange, 0.0096, &run, None, &[], T).unwrap();
    assert!(clusters.len() > 1, "{clusters:?}");
    for c in &clusters {
        assert!(c.span <= 0.06 + 1e-9, "{clusters:?}");
    }
}

/// A line touching a cylinder of radius 2 parts from it as `s^2 / 4`
/// along it. A crossing on the line found at the touch, by a section a
/// micron off the cylinder, may stand `sqrt(2 * 1e-6 * 2)` from the
/// true one either way; found `d` = 1e-3 along the line, the valley
/// reaches past the touch to `d + sqrt(d^2 + 4e-6)`.
/// A line crossing the cylinder at an angle, and an exact section,
/// spread nothing.
#[test]
fn a_crossing_where_the_edge_touches_the_other_surface_spreads_across_its_valley() {
    let cylinder: SurfaceGeometry = ogeom_geom::CylinderSurface::new(
        ogeom_math::Cylinder::new(Frame::WORLD, 2.0, T).unwrap(),
        (-5.0, 5.0),
    )
    .unwrap()
    .into();
    let line = |x: f64| {
        Curve::Line(
            ogeom_geom::LineCurve::segment(Point::new(x, -3.0, 0.0), Point::new(x, 3.0, 0.0), T)
                .unwrap(),
        )
    };
    let touching = line(2.0);
    let spread_at = |curve: &Curve, y: f64, honesty: f64| {
        let t = y + 3.0;
        let at = curve.point_at(t, T).unwrap();
        touch_spread((curve, t), &cylinder, at, honesty, T)
    };
    let at_touch = spread_at(&touching, 0.0, 1e-6);
    assert!((at_touch - 2e-3).abs() < 1e-9, "{at_touch}");
    let d = 1e-3;
    let beside = spread_at(&touching, d, 1e-6);
    let want = d + d.mul_add(d, 4e-6).sqrt();
    assert!((beside - want).abs() < 1e-6, "{beside} against {want}");
    assert_eq!(spread_at(&touching, 0.0, 0.0), 0.0);
    assert_eq!(spread_at(&line(1.0), 3.0_f64.sqrt(), 1e-6), 0.0);
}

/// A crossing of a loose edge spreads along it as its tolerance over
/// the sine of its lean out of the crossing face, held to ten
/// tolerances; a tight edge's does not; and where the face's normal
/// cannot be read, at a cone's apex, the spread is none.
#[test]
fn a_crossing_spreads_along_a_loose_edge_by_its_lean() {
    let mut model = Model::new();
    let solid = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let mut gathered = gather(&model, &solid, false, None, T).unwrap();
    let e = gathered
        .faces
        .iter_mut()
        .flat_map(|f| f.edges.iter_mut())
        .find(|e| matches!(&*e.curve, Curve::Line(_)))
        .unwrap();
    let t = f64::midpoint(e.crange.0, e.crange.1);
    let at = e.curve.point_at(t, T).unwrap();
    let along = e.curve.d1_at(t, T).unwrap();
    let along = along / along.magnitude();
    let square = if along.z.abs() < 0.9 {
        along.cross(ogeom_math::Vector::Z)
    } else {
        along.cross(ogeom_math::Vector::X)
    };
    let square = square / square.magnitude();
    let plane_leaning = |lean: f64| -> SurfaceGeometry {
        let normal = square * lean.cos() + along * lean.sin();
        let frame = Frame::new(
            at,
            Direction::new(normal, T).unwrap(),
            Direction::new(normal.cross(square).cross(normal), T)
                .or_else(|_| Direction::new(along, T))
                .unwrap(),
            T,
        )
        .unwrap();
        ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(frame)).into()
    };
    assert_eq!(lean_spread(e, t, &plane_leaning(0.3), at, T), 0.0);
    e.tolerance = 0.01;
    for lean in [0.3_f64, 1.2] {
        let spread = lean_spread(e, t, &plane_leaning(lean), at, T);
        assert!(
            (spread - 0.01 / lean.sin()).abs() < 1e-9,
            "{lean}: {spread}"
        );
    }
    let spread = lean_spread(e, t, &plane_leaning(0.01), at, T);
    assert!((spread - 0.1).abs() < 1e-12, "{spread}");
    let apex = ogeom_geom::ConeSurface::new(
        ogeom_math::Cone::new(
            Frame::new(
                at,
                Direction::new(square, T).unwrap(),
                Direction::new(along, T).unwrap(),
                T,
            )
            .unwrap(),
            0.0,
            0.5,
            T,
        )
        .unwrap(),
        (0.0, 5.0),
    )
    .unwrap();
    let spread = lean_spread(e, t, &SurfaceGeometry::Cone(apex), at, T);
    assert!(spread.is_finite() && spread >= 0.0, "{spread}");
}

/// A point on a face's surface inside its trim lies within it; one
/// beyond its boundary lies within it only to the margin asked.
#[test]
fn a_point_beside_a_face_lies_within_it_only_to_the_margin() {
    let mut model = Model::new();
    let solid = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let gathered = gather(&model, &solid, false, None, T).unwrap();
    let top = gathered
        .faces
        .iter()
        .find(|f| f.bound.low().unwrap().z > 9.9)
        .unwrap();
    assert!(within_face(top, Point::new(10.0, 10.0, 10.0), 1e-7, T).unwrap());
    let beside = Point::new(10.0, -0.05, 10.0);
    assert!(!within_face(top, beside, 1e-4, T).unwrap());
    assert!(within_face(top, beside, 0.06, T).unwrap());
}

/// A circle's count from its sag in closed form stands its chords within
/// the strand sag, short of the cap both counts share, and never above
/// the count measuring the chords reaches.
#[test]
fn a_circle_is_sampled_to_the_strand_sag_in_closed_form() {
    for radius in [0.01, 1.0, 40.0] {
        let circle = ogeom_math::Circle2::new(ogeom_math::Frame2::WORLD, radius, T).unwrap();
        let pcurve = PlanarCurve::Circle(ogeom_geom::Circle2d::new(circle));
        for (lo, hi) in [(0.0, 2.0 * PI), (1.0, 2.5)] {
            let count = circle_sag_count(&pcurve, lo, hi, 8).unwrap();
            let step = (hi - lo) / f64::from(u32::try_from(count).unwrap());
            assert!(
                count == 512 || radius * (1.0 - (step / 2.0).cos()) <= STRAND_SAG * (1.0 + 1e-9)
            );
            let (measured, _) = sag_count(&|t| pcurve.point_at(t, T), lo, hi, 8).unwrap();
            assert!(count <= measured, "{count} against {measured}");
        }
    }
}

/// A narrow wedge whose last edge stops short of the column it meets,
/// once by more than the weld and once by less than the snap: closed,
/// a point level with the opening and outside the wedge reads outside,
/// and one inside still reads inside.
#[test]
fn an_outline_closes_across_the_openings_its_ends_leave() {
    let top = PI / 2.0;
    for short in [4e-4, 2.6e-8] {
        let mut lines = vec![
            vec![Point2::new(0.0, 0.0), Point2::new(0.06, 0.0)],
            vec![Point2::new(0.06, 0.0), Point2::new(0.06, top)],
            vec![Point2::new(0.06, top - short), Point2::new(0.0, 0.0)],
        ];
        bridge_outline_gaps(&mut lines);
        let lines: Vec<&[Point2]> = lines.iter().map(Vec::as_slice).collect();
        // Level with the opening, a turn round from the wedge, so its
        // ray runs across the wedge.
        let level = Point2::new(-5.298, top - short / 2.0);
        assert!(!arrange::inside_many(&lines, level), "{short}");
        assert!(
            arrange::inside_many(&lines, Point2::new(0.05, 0.5)),
            "{short}"
        );
    }
}

fn frame_at(origin: Point) -> Frame {
    Frame::new(origin, Direction::Z, Direction::X, T).unwrap()
}

fn boxes(model: &mut Model) -> (Shape, Shape) {
    // Overlapping in the corner cube [1,2]^3.
    let a = make_box(model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let b = make_box(
        model,
        frame_at(Point::new(1.0, 1.0, 1.0)),
        (2.0, 2.0, 2.0),
        T,
    )
    .unwrap();
    (a.shape, b.shape)
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    let fine = Deflection {
        chord: 1e-3,
        ..Deflection::default()
    };
    volume_properties(model, shape, fine, T).unwrap().mass
}

fn assert_valid(model: &Model, shape: &Shape) {
    let diagnosis = check(model, shape, T).unwrap();
    assert!(
        diagnosis.is_valid(),
        "the result fails validity: {:?}",
        diagnosis.problems
    );
}

#[test]
fn fuse_of_overlapping_boxes_has_the_inclusion_exclusion_volume() {
    let mut model = Model::new();
    let (a, b) = boxes(&mut model);
    let fused = fuse(&mut model, &a, &b, T).unwrap();
    assert_valid(&model, &fused.shape);
    assert!((volume(&model, &fused.shape) - 15.0).abs() < 1e-9);
    assert_eq!(
        fused.history.modified(&a),
        std::slice::from_ref(&fused.shape)
    );
}

#[test]
fn common_of_overlapping_boxes_is_the_overlap_cube() {
    let mut model = Model::new();
    let (a, b) = boxes(&mut model);
    let result = common(&mut model, &a, &b, T).unwrap();
    assert_valid(&model, &result.shape);
    assert!((volume(&model, &result.shape) - 1.0).abs() < 1e-9);
}

#[test]
fn cut_removes_the_overlap_from_the_first_argument() {
    let mut model = Model::new();
    let (a, b) = boxes(&mut model);
    let result = cut(&mut model, &a, &b, T).unwrap();
    assert_valid(&model, &result.shape);
    assert!((volume(&model, &result.shape) - 7.0).abs() < 1e-9);
}

#[test]
fn section_of_overlapping_boxes_is_the_six_segment_seam() {
    let mut model = Model::new();
    let (a, b) = boxes(&mut model);
    let result = section(&mut model, &a, &b, T).unwrap();
    let edges = explore(&model, &result.shape, Filter::OfType(ShapeType::Edge)).unwrap();
    assert_eq!(edges.len(), 6);
    // One closed loop: six edges end to end on six shared vertices.
    let vertices = ogeom_topo::explore_unique(&model, &result.shape, ShapeType::Vertex).unwrap();
    assert_eq!(vertices.len(), 6);
}

#[test]
fn fuse_of_disjoint_boxes_is_a_compound_of_both() {
    let mut model = Model::new();
    let a = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let b = make_box(
        &mut model,
        frame_at(Point::new(5.0, 0.0, 0.0)),
        (2.0, 2.0, 2.0),
        T,
    )
    .unwrap();
    let fused = fuse(&mut model, &a.shape, &b.shape, T).unwrap();
    assert_eq!(model.kind_of(&fused.shape).unwrap(), ShapeType::Compound);
    assert!((volume(&model, &fused.shape) - 16.0).abs() < 1e-9);
}

#[test]
fn cutting_a_through_post_leaves_a_slab_with_a_hole() {
    let mut model = Model::new();
    let slab = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
    let post = make_box(
        &mut model,
        frame_at(Point::new(1.5, 1.5, -1.0)),
        (1.0, 1.0, 3.0),
        T,
    )
    .unwrap();
    let result = cut(&mut model, &slab.shape, &post.shape, T).unwrap();
    assert_valid(&model, &result.shape);
    assert!((volume(&model, &result.shape) - 15.0).abs() < 1e-9);
}

#[test]
fn common_with_a_contained_box_is_that_box() {
    let mut model = Model::new();
    let outer = make_box(&mut model, Frame::WORLD, (6.0, 6.0, 6.0), T).unwrap();
    let inner = make_box(
        &mut model,
        frame_at(Point::new(2.0, 2.0, 2.0)),
        (2.0, 2.0, 2.0),
        T,
    )
    .unwrap();
    let result = common(&mut model, &outer.shape, &inner.shape, T).unwrap();
    assert_valid(&model, &result.shape);
    assert!((volume(&model, &result.shape) - 8.0).abs() < 1e-9);
}

#[test]
fn cutting_everything_away_leaves_an_empty_compound() {
    let mut model = Model::new();
    let outer = make_box(&mut model, Frame::WORLD, (6.0, 6.0, 6.0), T).unwrap();
    let inner = make_box(
        &mut model,
        frame_at(Point::new(2.0, 2.0, 2.0)),
        (2.0, 2.0, 2.0),
        T,
    )
    .unwrap();
    let result = cut(&mut model, &inner.shape, &outer.shape, T).unwrap();
    assert_eq!(model.kind_of(&result.shape).unwrap(), ShapeType::Compound);
    assert!(
        explore(&model, &result.shape, Filter::OfType(ShapeType::Face))
            .unwrap()
            .is_empty()
    );
}

#[test]
fn drilling_a_box_leaves_a_cylindrical_hole() {
    // Box minus a through-post cylinder: the box's top and bottom faces
    // come back with *circular* holes, the hole's wall is the cylinder's
    // own surface with its material side flipped, and the cylinder's
    // seam and both section circles all split and sew for the shell to
    // close.
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
    let drill = make_cylinder(
        &mut model,
        frame_at(Point::new(2.0, 2.0, -1.0)),
        0.5,
        3.0,
        T,
    )
    .unwrap();
    let result = cut(&mut model, &block.shape, &drill.shape, T).unwrap();
    assert_valid(&model, &result.shape);
    let exact = 16.0 - PI * 0.25;
    let got = volume(&model, &result.shape);
    assert!(
        (got - exact).abs() / exact < 2e-3,
        "volume {got} against {exact}"
    );
}

/// A long rod of radius 5 at right angles to a long drill of radius 56,
/// the axes 46 apart, runs in through the drill's wall and out again:
/// its side meets the drill in two loops 64 apart. Both cut the rod, and
/// what the drill holds of it is the integral of the rod's slices below
/// the drill's wall.
#[test]
fn a_rod_through_a_long_drill_s_wall_is_cut_on_both_sides() {
    let (big, small, offset) = (56.0, 5.0, 46.0);
    let mut model = Model::new();
    let rod = make_cylinder(
        &mut model,
        Frame::new(
            Point::new(5.0, -1000.0, offset),
            Direction::Y,
            Direction::X,
            T,
        )
        .unwrap(),
        small,
        2000.0,
        T,
    )
    .unwrap();
    let drill = make_cylinder(
        &mut model,
        Frame::new(Point::new(-900.0, 0.0, 0.0), Direction::X, Direction::Z, T).unwrap(),
        big,
        1800.0,
        T,
    )
    .unwrap();
    // The slice of the rod at `y` lies in the drill below z = h(y): a
    // disc's segment, integrated by Simpson's rule.
    let segment = |y: f64| -> f64 {
        let d = (big * big - y * y).max(0.0).sqrt() - offset;
        if d <= -small {
            0.0
        } else if d >= small {
            PI * small * small
        } else {
            small * small * (-d / small).acos() + d * (small * small - d * d).sqrt()
        }
    };
    let reach = (big * big - (offset - small).powi(2)).sqrt();
    let steps = 20_000_u32;
    let h = 2.0 * reach / f64::from(steps);
    let held = (0..=steps)
        .map(|i| {
            let w = if i == 0 || i == steps {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            w * segment(-reach + h * f64::from(i))
        })
        .sum::<f64>()
        * h
        / 3.0;
    let whole = PI * small * small * 2000.0;
    let kept = cut(&mut model, &rod.shape, &drill.shape, T).unwrap();
    assert_valid(&model, &kept.shape);
    let inside = common(&mut model, &rod.shape, &drill.shape, T).unwrap();
    assert_valid(&model, &inside.shape);
    let got = volume(&model, &inside.shape);
    assert!(
        (got - held).abs() / held < 2e-3,
        "common {got} against {held}"
    );
    let got = volume(&model, &kept.shape);
    let exact = whole - held;
    assert!(
        (got - exact).abs() / exact < 2e-3,
        "cut {got} against {exact}"
    );
}

#[test]
fn common_of_a_box_and_a_cylinder_is_the_post_inside_it() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
    let post = make_cylinder(
        &mut model,
        frame_at(Point::new(2.0, 2.0, -1.0)),
        0.5,
        3.0,
        T,
    )
    .unwrap();
    let result = common(&mut model, &block.shape, &post.shape, T).unwrap();
    assert_valid(&model, &result.shape);
    let exact = PI * 0.25;
    let got = volume(&model, &result.shape);
    assert!(
        (got - exact).abs() / exact < 2e-3,
        "volume {got} against {exact}"
    );
}

#[test]
fn fusing_a_post_onto_a_slab_adds_what_stands_proud() {
    let mut model = Model::new();
    let slab = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
    let post = make_cylinder(
        &mut model,
        frame_at(Point::new(2.0, 2.0, -1.0)),
        0.5,
        3.0,
        T,
    )
    .unwrap();
    let result = fuse(&mut model, &slab.shape, &post.shape, T).unwrap();
    assert_valid(&model, &result.shape);
    let exact = PI.mul_add(0.25 * 3.0, 16.0) - PI * 0.25;
    let got = volume(&model, &result.shape);
    assert!(
        (got - exact).abs() / exact < 2e-3,
        "volume {got} against {exact}"
    );
}

#[test]
fn crossed_cylinders_run_through_the_marched_sections() {
    // No closed form exists for cylinder/cylinder: the sections are
    // marched and fitted, and every stage downstream (paves against
    // seams on both charts, splitting, classification, rebuilding with
    // fitted pcurves, sewing, meshing) has to work within the fit's
    // stated budget. Unequal radii keep the tangential branch points
    // away. The volumes have no easy closed form, so the operations are
    // held to each other: fuse = A + B - common and cut = A - common are
    // identities whatever the shapes.
    let mut model = Model::new();
    let upright = make_cylinder(&mut model, Frame::WORLD, 1.0, 4.0, T).unwrap();
    let across_frame =
        Frame::new(Point::new(-2.0, 0.0, 2.0), Direction::X, Direction::Y, T).unwrap();
    let across = make_cylinder(&mut model, across_frame, 0.6, 4.0, T).unwrap();

    let both = fuse(&mut model, &upright.shape, &across.shape, T).unwrap();
    assert_valid(&model, &both.shape);
    let shared = common(&mut model, &upright.shape, &across.shape, T).unwrap();
    assert_valid(&model, &shared.shape);
    let pierced = cut(&mut model, &upright.shape, &across.shape, T).unwrap();
    assert_valid(&model, &pierced.shape);

    let va = volume(&model, &upright.shape);
    let vb = volume(&model, &across.shape);
    let vf = volume(&model, &both.shape);
    let vc = volume(&model, &shared.shape);
    let vx = volume(&model, &pierced.shape);

    assert!(vc > 0.0 && vc < vb, "the overlap is real and partial: {vc}");
    assert!(
        (vf - (va + vb - vc)).abs() / vf < 2e-3,
        "fuse {vf} against A + B - common {}",
        va + vb - vc
    );
    assert!(
        (vx - (va - vc)).abs() / vx < 2e-3,
        "cut {vx} against A - common {}",
        va - vc
    );
}

#[test]
fn stacked_boxes_fuse_into_one_solid_and_the_contact_vanishes() {
    // Same-domain contact with opposed materials: the shared rectangle is
    // interior to the union and no face of the result may carry it.
    let mut model = Model::new();
    let lower = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 1.0), T).unwrap();
    let upper = make_box(
        &mut model,
        frame_at(Point::new(0.0, 0.0, 1.0)),
        (2.0, 2.0, 1.0),
        T,
    )
    .unwrap();
    let fused = fuse(&mut model, &lower.shape, &upper.shape, T).unwrap();
    assert_valid(&model, &fused.shape);
    assert_eq!(model.kind_of(&fused.shape).unwrap(), ShapeType::Solid);
    assert!((volume(&model, &fused.shape) - 8.0).abs() < 1e-9);
}

#[test]
fn a_small_box_on_a_big_one_fuses_with_a_partial_contact() {
    // The big top face splits into the contact rectangle (swallowed)
    // and the surround, which stays and must sew to the small box's
    // walls along the contact's edges.
    let mut model = Model::new();
    let big = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
    let small = make_box(
        &mut model,
        frame_at(Point::new(1.0, 1.0, 1.0)),
        (2.0, 2.0, 1.0),
        T,
    )
    .unwrap();
    let fused = fuse(&mut model, &big.shape, &small.shape, T).unwrap();
    assert_valid(&model, &fused.shape);
    assert!((volume(&model, &fused.shape) - 20.0).abs() < 1e-9);

    // Cutting the same pair removes nothing but the measure-zero contact:
    // the big box survives whole, its top face's contact piece intact.
    let mut model = Model::new();
    let big = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 1.0), T).unwrap();
    let small = make_box(
        &mut model,
        frame_at(Point::new(1.0, 1.0, 1.0)),
        (2.0, 2.0, 1.0),
        T,
    )
    .unwrap();
    let result = cut(&mut model, &big.shape, &small.shape, T).unwrap();
    assert_valid(&model, &result.shape);
    assert!((volume(&model, &result.shape) - 16.0).abs() < 1e-9);
}

#[test]
fn flush_walls_fuse_cut_and_meet_with_aligned_contact() {
    // Overlapping boxes sharing flush walls: same-domain contact with
    // *aligned* materials. A = [0,2]^3, B = [1,3]x[0,2]x[0,2]: the y and
    // z walls of the overlap are coplanar with aligned outward normals.
    let mut model = Model::new();
    let a = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let b = make_box(
        &mut model,
        frame_at(Point::new(1.0, 0.0, 0.0)),
        (2.0, 2.0, 2.0),
        T,
    )
    .unwrap();
    let fused = fuse(&mut model, &a.shape, &b.shape, T).unwrap();
    assert_valid(&model, &fused.shape);
    assert!((volume(&model, &fused.shape) - 12.0).abs() < 1e-9);

    let mut model = Model::new();
    let a = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let b = make_box(
        &mut model,
        frame_at(Point::new(1.0, 0.0, 0.0)),
        (2.0, 2.0, 2.0),
        T,
    )
    .unwrap();
    let shared = common(&mut model, &a.shape, &b.shape, T).unwrap();
    assert_valid(&model, &shared.shape);
    assert!((volume(&model, &shared.shape) - 4.0).abs() < 1e-9);

    let mut model = Model::new();
    let a = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let b = make_box(
        &mut model,
        frame_at(Point::new(1.0, 0.0, 0.0)),
        (2.0, 2.0, 2.0),
        T,
    )
    .unwrap();
    let cut_result = cut(&mut model, &a.shape, &b.shape, T).unwrap();
    assert_valid(&model, &cut_result.shape);
    assert!((volume(&model, &cut_result.shape) - 4.0).abs() < 1e-9);
}

#[test]
fn every_source_face_is_accounted_for_in_the_history() {
    let mut model = Model::new();
    let (a, b) = boxes(&mut model);
    let result = cut(&mut model, &a, &b, T).unwrap();
    for face in explore(&model, &a, Filter::OfType(ShapeType::Face)).unwrap() {
        assert!(
            result.history.is_affected(&face),
            "a face of the first argument vanished from the history"
        );
    }
}

/// A plate 10 thick, `rows` by `rows` holes of radius 1.5 drilled
/// through it 10 apart and 10 from its sides.
fn drilled_plate(model: &mut Model, rows: u32) -> Shape {
    let size = 10.0 * f64::from(rows) + 10.0;
    let block = make_box(model, Frame::WORLD, (size, size, 10.0), T)
        .unwrap()
        .shape;
    let mut pins = Vec::new();
    for i in 0..rows {
        for j in 0..rows {
            let at = Point::new(10.0 + 10.0 * f64::from(i), 10.0 + 10.0 * f64::from(j), -5.0);
            pins.push(
                make_cylinder(model, frame_at(at), 1.5, 20.0, T)
                    .unwrap()
                    .shape,
            );
        }
    }
    let pins = model.add_compound(&pins).unwrap();
    cut(model, &block, &pins, T).unwrap().shape
}

/// The operation with every face of both solids taken into the fuse.
fn with_every_face(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    keep_aside: [bool; 2],
    keep: impl Fn(&FacePiece) -> Option<bool>,
) -> Built {
    let fused = general_fuse_classified(model, a, b, None, T).unwrap();
    let kept: Vec<(usize, bool)> = fused
        .pieces
        .iter()
        .enumerate()
        .filter_map(|(i, p)| keep(p).map(|flip| (i, flip)))
        .collect();
    assemble_result(model, &fused, &kept, keep_aside, a, b, T).unwrap()
}

fn cut_keeps(p: &FacePiece) -> Option<bool> {
    match (p.from_a, p.state) {
        (true, PieceState::Out | PieceState::OnOpposed) => Some(false),
        (false, PieceState::In) => Some(true),
        _ => None,
    }
}

/// Whether `face`'s box stands clear of each of `tool`'s lumps' by a
/// millimetre.
fn clear_of(model: &Model, face: &Shape, tool: &Shape) -> bool {
    let face = ogeom_algo::face_bounds(model, face).unwrap().expanded(1.0);
    explore_unique(model, tool, ShapeType::Solid)
        .unwrap()
        .iter()
        .all(|lump| !face.intersects(&ogeom_algo::tight_bounds(model, lump, T).unwrap()))
}

/// The result with faces set aside against the one with every face
/// taken: the same volume and face count, valid, and every face of `a`
/// clear of `b` passed through as itself.
fn same_as_with_every_face(model: &Model, a: &Shape, b: &Shape, local: &Built, whole: &Built) {
    assert_valid(model, &local.shape);
    assert_valid(model, &whole.shape);
    let (v_local, v_whole) = (volume(model, &local.shape), volume(model, &whole.shape));
    assert!(
        (v_local - v_whole).abs() <= 1e-9 * v_whole.abs(),
        "{v_local} against {v_whole}"
    );
    let count = |shape: &Shape| explore_unique(model, shape, ShapeType::Face).unwrap().len();
    assert_eq!(count(&local.shape), count(&whole.shape));
    let mut passed = 0;
    for face in explore_unique(model, a, ShapeType::Face).unwrap() {
        if clear_of(model, &face, b) {
            let copy = local.history.copy_of(&face);
            assert!(
                copy.is_some_and(|c| c.is_same(&face)),
                "a face clear of the tool was not passed through"
            );
            passed += 1;
        }
    }
    assert!(passed > 0);
}

/// A tool of one lump along two edges of a drilled plate, its box
/// spanning the plate while its faces reach only the edges: the hole
/// walls miss every face of it and stand outside it, so they are set
/// aside and passed through as they were, the result the one with every
/// face taken.
#[test]
fn a_one_lump_tool_spanning_a_plate_leaves_the_holes_it_misses_as_they_were() {
    let mut model = Model::new();
    let plate = drilled_plate(&mut model, 3);
    let along_y = make_box(
        &mut model,
        frame_at(Point::new(-1.0, -1.0, 8.0)),
        (5.0, 42.0, 3.0),
        T,
    )
    .unwrap()
    .shape;
    let along_x = make_box(
        &mut model,
        frame_at(Point::new(-1.0, -1.0, 9.0)),
        (42.0, 5.0, 2.0),
        T,
    )
    .unwrap()
    .shape;
    let corner = fuse(&mut model, &along_y, &along_x, T).unwrap().shape;
    assert_eq!(
        explore_unique(&model, &corner, ShapeType::Solid)
            .unwrap()
            .len(),
        1
    );
    let local = cut(&mut model, &plate, &corner, T).unwrap();
    let whole = with_every_face(&mut model, &plate, &corner, [true, false], cut_keeps);
    assert_valid(&model, &local.shape);
    let (v_local, v_whole) = (volume(&model, &local.shape), volume(&model, &whole.shape));
    assert!(
        (v_local - v_whole).abs() <= 1e-9 * v_whole,
        "{v_local} against {v_whole}"
    );
    let walls: Vec<Shape> = explore_unique(&model, &plate, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::Cylinder(_))
            )
        })
        .collect();
    assert_eq!(walls.len(), 9);
    for wall in &walls {
        let copy = local.history.copy_of(wall);
        assert!(
            copy.is_some_and(|c| c.is_same(wall)),
            "a hole wall clear of the tool was not passed through"
        );
    }
}

#[test]
fn a_hole_drilled_at_a_plates_corner_leaves_the_rest_as_it_was() {
    let mut model = Model::new();
    let plate = drilled_plate(&mut model, 3);
    let corner = make_cylinder(
        &mut model,
        frame_at(Point::new(5.0, 5.0, -5.0)),
        1.5,
        20.0,
        T,
    )
    .unwrap()
    .shape;
    assert!(
        aside::set_aside(&model, &plate, &corner, T)
            .unwrap()
            .is_some()
    );
    let local = cut(&mut model, &plate, &corner, T).unwrap();
    let whole = with_every_face(&mut model, &plate, &corner, [true, false], cut_keeps);
    same_as_with_every_face(&model, &plate, &corner, &local, &whole);
    let holed = 40.0 * 40.0 * 10.0 - 10.0 * PI * 1.5 * 1.5 * 10.0;
    assert!((volume(&model, &local.shape) - holed).abs() < 1e-6 * holed);
}

/// Two pins at opposite corners, cut as one compound whose box spans
/// the plate: the faces clear of both pins are set aside all the same.
#[test]
fn holes_at_opposite_corners_leave_the_faces_between_them_as_they_were() {
    let mut model = Model::new();
    let plate = drilled_plate(&mut model, 3);
    let mut pins = Vec::new();
    for at in [5.0, 35.0] {
        pins.push(
            make_cylinder(&mut model, frame_at(Point::new(at, at, -5.0)), 1.5, 20.0, T)
                .unwrap()
                .shape,
        );
    }
    let pins = model.add_compound(&pins).unwrap();
    let local = cut(&mut model, &plate, &pins, T).unwrap();
    let whole = with_every_face(&mut model, &plate, &pins, [true, false], cut_keeps);
    same_as_with_every_face(&model, &plate, &pins, &local, &whole);
    let holed = 40.0 * 40.0 * 10.0 - 11.0 * PI * 1.5 * 1.5 * 10.0;
    assert!((volume(&model, &local.shape) - holed).abs() < 1e-6 * holed);
}

/// Four bars along a box's top edges, fused into a frame whose top
/// stays split where the bars met: each edge between two of its pieces
/// lies on the box's top and is laid there once, and the cut leaves the
/// box less the frame.
#[test]
fn a_frame_fused_from_bars_cuts_the_rim_off_a_box() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let bars = [
        ((0.0, 0.0), (10.0, 1.0)),
        ((0.0, 9.0), (10.0, 1.0)),
        ((0.0, 0.0), (1.0, 10.0)),
        ((9.0, 0.0), (1.0, 10.0)),
    ];
    let mut frame: Option<Shape> = None;
    for ((x, y), (dx, dy)) in bars {
        let bar = make_box(
            &mut model,
            frame_at(Point::new(x, y, 9.0)),
            (dx, dy, 1.0),
            T,
        )
        .unwrap()
        .shape;
        frame = Some(match frame {
            None => bar,
            Some(held) => fuse(&mut model, &held, &bar, T).unwrap().shape,
        });
    }
    let frame = frame.unwrap();
    let rimmed = cut(&mut model, &block, &frame, T).unwrap();
    assert_valid(&model, &rimmed.shape);
    let want = 1000.0 - (100.0 - 64.0);
    let got = volume(&model, &rimmed.shape);
    assert!((got - want).abs() < 1e-9 * want, "{got} against {want}");
}

#[test]
fn a_notch_through_a_plates_edge_meets_the_faces_left_as_they_were() {
    // A box notching the corner through two sides, the top and the
    // bottom: the outer rings of all four are split, and the edges
    // they share with the faces set aside stay as they were.
    let mut model = Model::new();
    let plate = drilled_plate(&mut model, 3);
    let notch = make_box(
        &mut model,
        frame_at(Point::new(-1.0, -1.0, -1.0)),
        (5.0, 4.0, 12.0),
        T,
    )
    .unwrap()
    .shape;
    let local = cut(&mut model, &plate, &notch, T).unwrap();
    let whole = with_every_face(&mut model, &plate, &notch, [true, false], cut_keeps);
    same_as_with_every_face(&model, &plate, &notch, &local, &whole);
}

#[test]
fn a_boss_fused_and_kept_in_common_beside_faces_set_aside() {
    let mut model = Model::new();
    let plate = drilled_plate(&mut model, 3);
    let boss = make_cylinder(
        &mut model,
        frame_at(Point::new(5.0, 5.0, 5.0)),
        3.0,
        10.0,
        T,
    )
    .unwrap()
    .shape;
    let local = fuse(&mut model, &plate, &boss, T).unwrap();
    let whole = with_every_face(&mut model, &plate, &boss, [true, true], |p| {
        (p.state == PieceState::Out || (p.state == PieceState::OnAligned && !p.covered))
            .then_some(false)
    });
    same_as_with_every_face(&model, &plate, &boss, &local, &whole);

    let local = common(&mut model, &plate, &boss, T).unwrap();
    let whole = with_every_face(&mut model, &plate, &boss, [false, false], |p| {
        (p.state == PieceState::In || (p.state == PieceState::OnAligned && !p.covered))
            .then_some(false)
    });
    assert_valid(&model, &local.shape);
    let (v_local, v_whole) = (volume(&model, &local.shape), volume(&model, &whole.shape));
    assert!(
        (v_local - v_whole).abs() <= 1e-9 * v_whole,
        "{v_local} against {v_whole}"
    );
    for face in explore_unique(&model, &plate, ShapeType::Face).unwrap() {
        if clear_of(&model, &face, &boss) {
            assert!(local.history.is_deleted(&face));
        }
    }
}

/// The corner hole into a plate of 582 faces, with every face but the
/// top, the bottom and the tool's set aside, comes out the solid the fuse
/// of every face makes, and passes every face it leaves alone through as
/// itself.
#[test]
#[ignore = "heavy"]
fn a_corner_hole_in_a_large_plate_is_the_solid_every_face_makes() {
    let mut model = Model::new();
    let plate = drilled_plate(&mut model, 24);
    assert_eq!(
        explore_unique(&model, &plate, ShapeType::Face)
            .unwrap()
            .len(),
        582
    );
    let corner = make_cylinder(
        &mut model,
        frame_at(Point::new(5.0, 5.0, -5.0)),
        1.5,
        20.0,
        T,
    )
    .unwrap()
    .shape;
    let local = cut(&mut model, &plate, &corner, T).unwrap();
    let whole = with_every_face(&mut model, &plate, &corner, [true, false], cut_keeps);
    same_as_with_every_face(&model, &plate, &corner, &local, &whole);
    // The top and the bottom take the new hole; every other face is
    // untouched.
    let copied = explore_unique(&model, &plate, ShapeType::Face)
        .unwrap()
        .iter()
        .filter(|face| local.history.copy_of(face).is_some())
        .count();
    assert_eq!(copied, 580);
}

/// The face a primitive builder gave `role`.
fn face_with_role(model: &Model, solid: &Shape, role: ogeom_core::Role) -> Shape {
    explore_unique(model, solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            matches!(
                model.provenance_of(f),
                Some(
                    ogeom_core::Provenance::Primitive { role: r, .. }
                        | ogeom_core::Provenance::Derived { role: r, .. }
                ) if *r == role
            )
        })
        .unwrap()
}

/// The faces of `result` that `history` traces from `source`.
fn traced_into(model: &Model, result: &Shape, history: &History, source: &Shape) -> Vec<Shape> {
    explore_unique(model, result, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|r| history.trace(source).iter().any(|t| t.is_same(r)))
        .collect()
}

/// A cylinder's caps lying in a cube's faces: the common keeps the
/// cube's piece of each, and each cap traces to that face; the fuse
/// keeps one face per cap plane too, which both trace to. The
/// cylinder's side, half outside the cube, is not lost in the common.
#[test]
fn a_coincident_piece_dropped_for_the_others_traces_to_it() {
    use ogeom_algo::primitive::roles::{FACE_MAX_Z, FACE_MIN_Z};
    let mut model = Model::new();
    let cube = make_box(&mut model, Frame::WORLD, (100.0, 100.0, 100.0), T)
        .unwrap()
        .shape;
    let cylinder = make_cylinder(&mut model, Frame::WORLD, 50.0, 100.0, T)
        .unwrap()
        .shape;
    for keep_in in [true, false] {
        let built = if keep_in {
            common(&mut model, &cube, &cylinder, T).unwrap()
        } else {
            fuse(&mut model, &cube, &cylinder, T).unwrap()
        };
        let expected = if keep_in {
            PI * 50.0 * 50.0 * 100.0 / 4.0
        } else {
            1e6 + PI * 50.0 * 50.0 * 100.0 * 3.0 / 4.0
        };
        assert_valid(&model, &built.shape);
        let v = volume(&model, &built.shape);
        assert!(
            (v - expected).abs() <= 1e-3 * expected,
            "{v} against {expected}"
        );
        for role in [FACE_MIN_Z, FACE_MAX_Z] {
            let (side, cap) = (
                face_with_role(&model, &cube, role),
                face_with_role(&model, &cylinder, role),
            );
            assert!(!built.history.is_deleted(&cap));
            let from_cap = traced_into(&model, &built.shape, &built.history, &cap);
            let from_side = traced_into(&model, &built.shape, &built.history, &side);
            if keep_in {
                // The common keeps one face per cap plane, the one both
                // sources trace to and nothing else.
                assert_eq!(from_cap.len(), 1, "{role:?}: the cap");
                assert_eq!(from_side.len(), 1, "{role:?}: the cube's face");
            }
            let shared = from_cap
                .iter()
                .filter(|f| from_side.iter().any(|g| g.is_same(f)))
                .count();
            assert_eq!(shared, 1, "common {keep_in}");
        }
    }
}

/// A slot cut into a step already cut, its floor partly on the step's
/// floor: the cut keeps the step's floor there and the slot's floor
/// traces to it as well as to its own piece; the slot's ends, which lie
/// in the block's sides with the material on the same side, are cut
/// away and stay deleted.
#[test]
fn a_tool_face_on_a_floor_a_cut_left_traces_to_that_floor() {
    use ogeom_algo::primitive::roles::{FACE_MAX_Z, FACE_MIN_X, FACE_MIN_Z};
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (100.0, 100.0, 50.0), T)
        .unwrap()
        .shape;
    let step = make_box(
        &mut model,
        frame_at(Point::new(50.0, -10.0, 30.0)),
        (60.0, 120.0, 20.0),
        T,
    )
    .unwrap()
    .shape;
    let slot = make_box(
        &mut model,
        frame_at(Point::new(-10.0, 40.0, 30.0)),
        (120.0, 20.0, 30.0),
        T,
    )
    .unwrap()
    .shape;
    let (step_floor, slot_floor) = (
        face_with_role(&model, &step, FACE_MIN_Z),
        face_with_role(&model, &slot, FACE_MIN_Z),
    );
    let stepped = cut(&mut model, &block, &step, T).unwrap();
    let slotted = cut(&mut model, &stepped.shape, &slot, T).unwrap();
    assert_valid(&model, &slotted.shape);
    let history = stepped.history.then(&slotted.history);
    let v = volume(&model, &slotted.shape);
    assert!((v - 380_000.0).abs() <= 1e-6 * 380_000.0, "{v}");
    let from_slot = traced_into(&model, &slotted.shape, &history, &slot_floor);
    let from_step = traced_into(&model, &slotted.shape, &history, &step_floor);
    assert_eq!(from_slot.len(), 2);
    let shared = from_slot
        .iter()
        .filter(|f| from_step.iter().any(|g| g.is_same(f)))
        .count();
    assert_eq!(shared, 1);
    // Faces whose material is gone stay gone: the slot's top and its
    // end outside the block.
    for role in [FACE_MAX_Z, FACE_MIN_X] {
        assert!(history.is_deleted(&face_with_role(&model, &slot, role)));
    }
}

/// A cut followed by a move: each face of the cut's inputs that is not
/// deleted traces to faces of the moved result, and to nothing else.
#[test]
fn a_cut_then_a_move_traces_its_inputs_to_the_moved_result() {
    let mut model = Model::new();
    let (a, b) = boxes(&mut model);
    let cut_out = cut(&mut model, &a, &b, T).unwrap();
    let moved = ogeom_algo::transformed(
        &mut model,
        &cut_out.shape,
        ogeom_math::Transform::translation(ogeom_math::Vector::new(100.0, 0.0, 0.0)),
    )
    .unwrap();
    let history = cut_out.history.then(&moved.history);
    let result = explore_unique(&model, &moved.shape, ShapeType::Face).unwrap();
    let mut traced = 0;
    for face in explore_unique(&model, &a, ShapeType::Face)
        .unwrap()
        .into_iter()
        .chain(explore_unique(&model, &b, ShapeType::Face).unwrap())
    {
        if history.is_deleted(&face) {
            continue;
        }
        let images = history.trace(&face);
        assert!(!images.is_empty());
        for image in images {
            assert!(result.iter().any(|r| r.is_same(image)));
        }
        traced += 1;
    }
    assert_eq!(traced, 9);
}

/// A cut followed by a shear: each face, edge and vertex of the cut's
/// inputs that the cut carries into its result traces through the
/// composed history to as many sub-shapes of its kind in the sheared
/// result, and to nothing else.
#[test]
fn a_cut_then_a_shear_traces_its_inputs_to_the_sheared_result() {
    let mut model = Model::new();
    let (a, b) = boxes(&mut model);
    let cut_out = cut(&mut model, &a, &b, T).unwrap();
    let shear = ogeom_math::GeneralTransform::new(
        ogeom_math::Matrix3 {
            rows: [[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        },
        ogeom_math::Vector::ZERO,
    );
    let sheared =
        ogeom_algo::general_transformed_shape(&mut model, &cut_out.shape, &shear, T).unwrap();
    let history = cut_out.history.then(&sheared.history);
    for (kind, expected) in [
        (ShapeType::Face, 9),
        (ShapeType::Edge, 15),
        (ShapeType::Vertex, 8),
    ] {
        let between = explore_unique(&model, &cut_out.shape, kind).unwrap();
        let result = explore_unique(&model, &sheared.shape, kind).unwrap();
        let mut traced = 0;
        for part in explore_unique(&model, &a, kind)
            .unwrap()
            .into_iter()
            .chain(explore_unique(&model, &b, kind).unwrap())
        {
            let carried = cut_out.history.trace(&part);
            if carried.is_empty() || !carried.iter().all(|c| between.iter().any(|r| r.is_same(c))) {
                continue;
            }
            let images = history.trace(&part);
            assert_eq!(images.len(), carried.len(), "{kind:?}");
            assert!(history.is_affected(&part), "{kind:?}");
            for image in images {
                assert!(result.iter().any(|r| r.is_same(image)), "{kind:?}");
            }
            traced += 1;
        }
        assert_eq!(traced, expected, "{kind:?}");
    }
}

/// Each edge and vertex of `argument` traced through `history`: deleted,
/// or traced to shapes of its kind that `result` holds. The counts of
/// those deleted, unchanged (held by the result as they stand, or
/// copied) and modified, by kind.
fn sub_shapes_traced(
    model: &Model,
    argument: &Shape,
    result: &Shape,
    history: &History,
) -> [(usize, usize, usize); 2] {
    let mut counts = [(0, 0, 0); 2];
    for (slot, kind) in [ShapeType::Edge, ShapeType::Vertex].into_iter().enumerate() {
        let held = explore_unique(model, result, kind).unwrap();
        for part in explore_unique(model, argument, kind).unwrap() {
            if history.is_deleted(&part) {
                counts[slot].0 += 1;
                continue;
            }
            let images = history.trace(&part);
            assert!(!images.is_empty(), "{kind:?}");
            for image in images {
                assert!(
                    held.iter().any(|h| h.is_same(image)),
                    "a {kind:?} traces to a shape the result does not hold"
                );
            }
            if history.copy_of(&part).is_some() || !history.is_affected(&part) {
                counts[slot].1 += 1;
            } else {
                counts[slot].2 += 1;
            }
        }
    }
    counts
}

fn vertex_point(model: &Model, vertex: &Shape) -> Point {
    let data = model.node(vertex).unwrap().data().as_vertex().unwrap();
    vertex.transform(model.datums()).unwrap().apply(data.point)
}

/// A drill through a block whose top line, its seam, stands a hair
/// above the block's top face: the face is cut along two lines either
/// side of the seam, closer than the weld, which sew into one line the
/// drill's face meets on both sides. That line is the drill's seam, and
/// the result is valid with the volume the drill exactly tangent to the
/// top cuts.
#[test]
fn a_drill_whose_seam_grazes_a_face_keeps_it_a_seam() {
    let chord = Deflection::with_chord(0.01).unwrap();
    let volume_at = |rise: f64| {
        let mut model = Model::new();
        let block = make_box(
            &mut model,
            frame_at(Point::new(-20.0, 0.0, 0.0)),
            (40.0, 30.0, 10.0),
            T,
        )
        .unwrap()
        .shape;
        let up = ogeom_math::Direction::Y;
        let seam = ogeom_math::Direction::Z;
        let drill = make_cylinder(
            &mut model,
            Frame::new(Point::new(0.0, -5.0, 7.0 + rise), up, seam, T).unwrap(),
            3.0,
            40.0,
            T,
        )
        .unwrap()
        .shape;
        let out = cut(&mut model, &block, &drill, T).unwrap();
        let d = check(&model, &out.shape, T).unwrap();
        assert!(d.is_valid(), "rise {rise}: {d}");
        volume_properties(&model, &out.shape, chord, T)
            .unwrap()
            .mass
    };
    let tangent = volume_at(0.0);
    for rise in [1e-9, 1e-8, 1e-7, 2e-7] {
        let v = volume_at(rise);
        assert!(
            (v - tangent).abs() < 1e-6 * tangent,
            "rise {rise}: {v} against {tangent}"
        );
    }
}

/// A drill tangent to a wall along the edge where a concave round
/// leaves it, the drill flatter than the round: the drill touches the
/// wall from the empty side and cuts into the material behind the
/// round, so its wall lines the cut from the tangent line to where it
/// leaves the solid. The material is a slab with a notch whose inner
/// corner is rounded, and the volume removed is the drill's disc over
/// the material's section, in closed form.
#[test]
fn a_drill_tangent_to_a_wall_along_a_rounds_edge_cuts_the_sliver_behind_it() {
    let mut model = Model::new();
    let slab = make_box(&mut model, Frame::WORLD, (40.0, 40.0, 10.0), T)
        .unwrap()
        .shape;
    let notch = make_box(
        &mut model,
        frame_at(Point::new(20.0, 20.0, -1.0)),
        (21.0, 21.0, 12.0),
        T,
    )
    .unwrap()
    .shape;
    let notched = cut(&mut model, &slab, &notch, T).unwrap().shape;
    // The round: the corner square less a cylinder of radius 2 about
    // (22, 22), tangent to the walls y = 20 and x = 20 at (22, 20) and
    // (20, 22).
    let square = make_box(
        &mut model,
        frame_at(Point::new(20.0, 20.0, 0.0)),
        (2.0, 2.0, 10.0),
        T,
    )
    .unwrap()
    .shape;
    let bore = make_cylinder(
        &mut model,
        frame_at(Point::new(22.0, 22.0, -1.0)),
        2.0,
        12.0,
        T,
    )
    .unwrap()
    .shape;
    let round = cut(&mut model, &square, &bore, T).unwrap().shape;
    let part = fuse(&mut model, &notched, &round, T).unwrap().shape;
    let chord = Deflection::with_chord(0.01).unwrap();
    let before = volume_properties(&model, &part, chord, T).unwrap().mass;
    for seam in [
        Direction::X,
        Direction::new(ogeom_math::Vector::new(-1.0, 0.0, 0.0), T).unwrap(),
    ] {
        let drill = make_cylinder(
            &mut model,
            Frame::new(Point::new(22.0, 50.0, -1.0), Direction::Z, seam, T).unwrap(),
            30.0,
            12.0,
            T,
        )
        .unwrap()
        .shape;
        let out = cut(&mut model, &part, &drill, T).unwrap();
        check(&model, &out.shape, T).unwrap();
        // The integral of sqrt(r^2 - u^2).
        let seg = |u: f64, r: f64| 0.5 * (u * (r * r - u * u).sqrt() + r * r * (u / r).asin());
        // Over the block x in [0, 20]: between y = 40 and the drill's
        // lower arc, y = 50 - sqrt(900 - (x - 22)^2).
        let block = (seg(-2.0, 30.0) - seg(-22.0, 30.0)) - 10.0 * 20.0;
        // Over the corner square x in [20, 22]: the square less its
        // quarter disc, less what lies below the drill's arc.
        let corner =
            (4.0 - (seg(0.0, 2.0) - seg(-2.0, 2.0))) - (60.0 - (seg(0.0, 30.0) - seg(-2.0, 30.0)));
        let expected = before - 10.0 * (block + corner);
        let after = volume_properties(&model, &out.shape, chord, T)
            .unwrap()
            .mass;
        assert!(
            (after - expected).abs() < 1e-6 * expected,
            "seam {seam:?}: {after} against {expected}"
        );
    }
}

/// A box cut by one overlapping its corner: the target's three edges
/// through the corner are shortened and the corner deleted; the tool's
/// three edges through its corner inside the target line the cavity,
/// cut short, and the rest of the tool is deleted.
#[test]
fn a_cut_records_every_edge_and_vertex_of_both_boxes() {
    let mut model = Model::new();
    let (a, b) = boxes(&mut model);
    let out = cut(&mut model, &a, &b, T).unwrap();
    let h = &out.history;
    // Edges (deleted, copied, modified), then vertices.
    assert_eq!(
        sub_shapes_traced(&model, &a, &out.shape, h),
        [(0, 9, 3), (1, 7, 0)]
    );
    assert_eq!(
        sub_shapes_traced(&model, &b, &out.shape, h),
        [(9, 0, 3), (7, 0, 1)]
    );
    for vertex in explore_unique(&model, &a, ShapeType::Vertex).unwrap() {
        let at = vertex_point(&model, &vertex);
        assert_eq!(
            h.is_deleted(&vertex),
            at.distance(Point::new(2.0, 2.0, 2.0)) < 1e-9
        );
    }
    // Each shortened edge is one edge of length one.
    for edge in explore_unique(&model, &a, ShapeType::Edge).unwrap() {
        if h.copy_of(&edge).is_some() || !h.is_affected(&edge) {
            continue;
        }
        let [image] = h.modified(&edge) else {
            panic!("a shortened edge is one edge");
        };
        let ends: Vec<Point> = explore_unique(&model, image, ShapeType::Vertex)
            .unwrap()
            .iter()
            .map(|v| vertex_point(&model, v))
            .collect();
        assert_eq!(ends.len(), 2);
        assert!((ends[0].distance(ends[1]) - 1.0).abs() < 1e-9);
    }
}

/// A plate drilled through: its edges and vertices come through as
/// copies, the drill's side seam is cut down to the plate's thickness,
/// and its rims and their vertices, outside the plate, are deleted.
#[test]
fn a_drill_records_every_edge_and_vertex_of_plate_and_drill() {
    let mut model = Model::new();
    let plate = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let drill = make_cylinder(
        &mut model,
        frame_at(Point::new(10.0, 10.0, -1.0)),
        3.0,
        12.0,
        T,
    )
    .unwrap()
    .shape;
    let out = cut(&mut model, &plate, &drill, T).unwrap();
    let h = &out.history;
    assert_eq!(
        sub_shapes_traced(&model, &plate, &out.shape, h),
        [(0, 12, 0), (0, 8, 0)]
    );
    assert_eq!(
        sub_shapes_traced(&model, &drill, &out.shape, h),
        [(2, 0, 1), (2, 0, 0)]
    );
}

/// A cut composed with a move: what the cut records of each edge and
/// vertex is carried to the moved result, and a deletion stays one.
#[test]
fn a_cut_then_a_move_traces_every_edge_and_vertex() {
    let mut model = Model::new();
    let (a, b) = boxes(&mut model);
    let out = cut(&mut model, &a, &b, T).unwrap();
    let moved = ogeom_algo::transformed(
        &mut model,
        &out.shape,
        ogeom_math::Transform::translation(ogeom_math::Vector::new(5.0, 0.0, 0.0)),
    )
    .unwrap();
    let composed = out.history.then(&moved.history);
    // A move modifies what it moves: the cut's copies are modified
    // in the composed history, and the cut's deletions stay deleted.
    let carried = |counts: [(usize, usize, usize); 2]| counts.map(|(d, c, m)| (d, c + m));
    for argument in [&a, &b] {
        let composed_counts = sub_shapes_traced(&model, argument, &moved.shape, &composed);
        assert!(composed_counts.iter().all(|&(_, copied, _)| copied == 0));
        assert_eq!(
            carried(composed_counts),
            carried(sub_shapes_traced(
                &model,
                argument,
                &out.shape,
                &out.history
            )),
        );
        for kind in [ShapeType::Edge, ShapeType::Vertex] {
            for part in explore_unique(&model, argument, kind).unwrap() {
                assert_eq!(composed.trace(&part).len(), out.history.trace(&part).len());
            }
        }
    }
}
