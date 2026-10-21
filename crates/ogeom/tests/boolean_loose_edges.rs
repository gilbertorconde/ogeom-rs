//! A face with one loose edge welds what is near that edge, and a section
//! crossing another, tight, edge of the face near its vertex keeps the
//! stretch it cuts there, as the face across that edge does.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::{Tolerance, Tolerances};
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The loose edge's tolerance, as a converted mesh's threaded edge states.
const LOOSE: f64 = 0.017;

/// Measured at a chord below the sliver's width: a coarser one draws the
/// sliver's thin faces across it and overcounts it.
fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::with_chord(1e-5).unwrap(), T)
        .unwrap()
        .mass
}

fn point_of(model: &Model, vertex: &Shape) -> Point {
    match model.node(vertex).unwrap().data() {
        NodeData::Vertex(data) => data.point,
        _ => unreachable!("a vertex"),
    }
}

/// A 20 x 20 x 10 box whose front face's top edge (y = 0, z = 10) states
/// `LOOSE`, its two vertices with it or left as built.
fn box_with_a_loose_top_edge(model: &mut Model, with_vertices: bool) -> Shape {
    let solid = ogeom::algo::make_box(model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let mut widened = 0;
    for edge in explore_unique(model, &solid, ShapeType::Edge).unwrap() {
        let ends = explore_unique(model, &edge, ShapeType::Vertex).unwrap();
        if !ends.iter().all(|v| {
            let p = point_of(model, v);
            p.y.abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9
        }) {
            continue;
        }
        widened += 1;
        if with_vertices {
            model.widen(&edge, Tolerance::new(LOOSE).unwrap()).unwrap();
        } else if let Some(node) = model.node_mut(&edge)
            && let NodeData::Edge(data) = node.data_mut()
        {
            data.tolerance = data.tolerance.widen_to(LOOSE);
        }
    }
    assert_eq!(widened, 1);
    solid
}

/// The area of the disc of `radius` about `(cx, cy)` in the quadrant
/// x >= 0, y >= 0, by the midpoint rule over its vertical chords.
fn disc_in_quadrant(cx: f64, cy: f64, radius: f64) -> f64 {
    let (lo, hi) = ((cx - radius).max(0.0), cx + radius);
    let steps = 100_000;
    let h = (hi - lo) / f64::from(steps);
    (0..steps)
        .map(|i| {
            let x = (f64::from(i) + 0.5).mul_add(h, lo);
            let half = (x - cx).mul_add(-(x - cx), radius * radius).max(0.0).sqrt();
            ((cy + half) - (cy - half).max(0.0)).max(0.0) * h
        })
        .sum()
}

/// A vertical drill of radius sqrt(1.25) about (1 + d, 0.5) crosses the
/// front face's tight bottom edge `d` from the corner vertex (0, 0, 0) and
/// the loose top edge `d` from the vertex above it, and leaves a sliver of
/// the box in the corner between them. The loose edge welds the crossings
/// near it, and the crossing of the tight edge stays where it is on every
/// face, so the cut and the common check valid. Their volumes are the box
/// less the drill's share of it, and that share, to within the sliver and
/// what the weld at the top may move: a vertex moved by up to the loose
/// edge's doubt at the top of a face the box's height tall. Together they
/// make the block.
fn drill_beside_the_corner(d: f64, with_vertices: bool) {
    let mut model = Model::with_tolerances(T);
    let solid = box_with_a_loose_top_edge(&mut model, with_vertices);
    let radius = 1.25_f64.sqrt();
    let at = Frame::new(
        Point::new(1.0 + d, 0.5, -5.0),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, at, radius, 20.0, T)
        .unwrap()
        .shape;
    let mut volumes = Vec::new();
    for name in ["cut", "common"] {
        let made = if name == "cut" {
            ogeom::boolean::cut(&mut model, &solid, &drill, T)
        } else {
            ogeom::boolean::common(&mut model, &solid, &drill, T)
        }
        .unwrap_or_else(|e| panic!("{name} at d = {d} (vertices {with_vertices}): {e}"))
        .shape;
        let diagnosis = ogeom::algo::check(&model, &made, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "{name} at d = {d} (vertices {with_vertices}): {diagnosis}"
        );
        volumes.push(volume(&model, &made));
    }
    let share = disc_in_quadrant(1.0 + d, 0.5, radius) * 10.0;
    // The sliver lies within d of the front face and within the drill's
    // crossing of the side face, y = 0.5 - sqrt(1.25 - (1 + d)^2), of the
    // other, the box's whole height.
    let across = 0.5 - (1.0 + d).mul_add(-(1.0 + d), 1.25).sqrt();
    let doubt = 2.0 * LOOSE;
    let allowed = d.mul_add(across * 10.0, doubt * doubt * 10.0 / 2.0);
    for (name, got, want) in [
        ("cut", volumes[0], 4000.0 - share),
        ("common", volumes[1], share),
    ] {
        assert!(
            (got - want).abs() <= allowed,
            "{name} at d = {d} (vertices {with_vertices}): {got} against {want}"
        );
    }
    let total = volumes[0] + volumes[1];
    assert!(
        (total - 4000.0).abs() < 4000.0 * 1e-6,
        "cut and common at d = {d} (vertices {with_vertices}): {total}"
    );
}

/// The issue's crossings, and others from a hundredth of a micron off the
/// vertex to past the loose edge's doubt.
#[test]
fn a_drill_beside_the_vertex_of_a_tight_edge_cuts_valid_by_a_loose_edge() {
    for with_vertices in [false, true] {
        for d in [1e-5, 0.001, 0.0065, 0.015, 0.02, 0.05] {
            drill_beside_the_corner(d, with_vertices);
        }
    }
}

/// A 20 x 20 x 10 box whose front top edge (y = 0, z = 10) stops `short`
/// before the corner (20, 0, 10) where the right top and front right edges
/// meet it, and whose vertex there stands midway between the two ends with
/// `radius`: the curves a file closes on one vertex stop either side of it,
/// each within its radius, as a fitted export's may.
fn box_whose_top_edge_stops_short(model: &mut Model, short: f64, radius: f64) -> Shape {
    use ogeom::geom::Curve3d as _;
    use ogeom::topo::EdgeRepr;
    let solid = ogeom::algo::make_box(model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let corner = Point::new(20.0, 0.0, 10.0);
    let mut moved = 0;
    for edge in explore_unique(model, &solid, ShapeType::Edge).unwrap() {
        let ends = explore_unique(model, &edge, ShapeType::Vertex).unwrap();
        if !ends.iter().all(|v| {
            let p = point_of(model, v);
            p.y.abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9
        }) {
            continue;
        }
        let NodeData::Edge(data) = model.node(&edge).unwrap().data() else {
            unreachable!("an edge")
        };
        let Some(&EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            unreachable!("a box edge has a curve")
        };
        let line = model.geometry().curve(curve).unwrap();
        let (from, to) = (
            line.point_at(range.0, T).unwrap(),
            line.point_at(range.1, T).unwrap(),
        );
        let by = (range.1 - range.0) * short / from.distance(to);
        let at_end = to.distance(corner) < 1e-9;
        assert!(at_end || from.distance(corner) < 1e-9);
        let shorten = |r: (f64, f64)| {
            if at_end {
                (r.0, r.1 - by)
            } else {
                (r.0 + by, r.1)
            }
        };
        let Some(node) = model.node_mut(&edge) else {
            unreachable!("the edge")
        };
        let NodeData::Edge(data) = node.data_mut() else {
            unreachable!("an edge")
        };
        for repr in &mut data.representations {
            match repr {
                EdgeRepr::Curve3d { range, .. }
                | EdgeRepr::PCurve { range, .. }
                | EdgeRepr::Seam { range, .. } => *range = shorten(*range),
                _ => {}
            }
        }
        moved += 1;
        for v in ends {
            if point_of(model, &v).distance(corner) < 1e-9
                && let Some(node) = model.node_mut(&v)
                && let NodeData::Vertex(data) = node.data_mut()
            {
                data.point = Point::new(20.0 - short / 2.0, 0.0, 10.0);
                data.tolerance = Tolerance::new(radius).unwrap();
            }
        }
    }
    assert_eq!(moved, 1);
    solid
}

/// The top edge stops 0.3 µm short of a corner whose vertex, midway,
/// holds 0.2 µm: each end lies within the vertex, but they lie further
/// apart than either end's own radius reaches. A pocket cut from the top
/// face far from that corner sections the face, and its boundary closes
/// at the vertex as the topology says. Cut and common hold the box less the
/// pocket and the pocket.
#[test]
fn a_joint_whose_curves_stop_either_side_of_its_vertex_closes_within_the_vertex() {
    let mut model = Model::with_tolerances(T);
    let solid = box_whose_top_edge_stops_short(&mut model, 3e-4, 2e-4);
    let at = Frame::new(Point::new(2.0, 2.0, 8.0), Direction::Z, Direction::X, T).unwrap();
    let pocket = ogeom::algo::make_box(&mut model, at, (4.0, 4.0, 4.0), T)
        .unwrap()
        .shape;
    for (name, want) in [("cut", 4000.0 - 32.0), ("common", 32.0)] {
        let made = if name == "cut" {
            ogeom::boolean::cut(&mut model, &solid, &pocket, T)
        } else {
            ogeom::boolean::common(&mut model, &solid, &pocket, T)
        }
        .unwrap_or_else(|e| panic!("{name}: {e}"))
        .shape;
        let diagnosis = ogeom::algo::check(&model, &made, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
        let got = volume(&model, &made);
        assert!(
            (got - want).abs() < 1e-6 * 4000.0,
            "{name}: {got} against {want}"
        );
    }
}

/// The solid `solid` built anew with the edge `dropped` left out and its
/// two vertices made one at their midpoint, holding both within its radius.
fn with_edge_collapsed(model: &mut Model, solid: &Shape, dropped: &Shape) -> Shape {
    use ogeom::core::FastMap;
    use ogeom::topo::{TShapeId, VertexData};
    let ends = explore_unique(model, dropped, ShapeType::Vertex).unwrap();
    let (a, b) = (point_of(model, &ends[0]), point_of(model, &ends[1]));
    let joined = model.add_vertex(
        VertexData::with_tolerance(Point::midpoint(a, b), a.distance(b) / 2.0 + T.confusion())
            .unwrap(),
    );
    let mut made: FastMap<TShapeId, Shape> = FastMap::default();
    for end in &ends {
        made.insert(end.node(), joined.clone());
    }
    fn rebuild(
        model: &mut Model,
        shape: &Shape,
        dropped: TShapeId,
        made: &mut FastMap<TShapeId, Shape>,
    ) -> Shape {
        if let Some(done) = made.get(&shape.node()) {
            return done.clone();
        }
        let node = model.node(shape).unwrap();
        let data = node.data().clone();
        let children: Vec<Shape> = node.children().to_vec();
        let mut kept = Vec::new();
        for child in children.iter().filter(|c| c.node() != dropped) {
            let bare = child
                .located(ogeom::topo::Location::identity())
                .oriented(ogeom::topo::Orientation::Forward);
            let new = rebuild(model, &bare, dropped, made);
            kept.push(
                new.located(child.location().clone())
                    .oriented(child.orientation()),
            );
        }
        let new = match (data, model.kind_of(shape).unwrap()) {
            (NodeData::Vertex(v), _) => model.add_vertex(v),
            (NodeData::Edge(e), _) => model.add_edge(*e, &kept).unwrap(),
            (NodeData::Face(f), _) => model.add_face(*f, &kept).unwrap(),
            (_, ShapeType::Wire) => model.add_wire(&kept).unwrap(),
            (_, ShapeType::Shell) => model.add_shell(&kept).unwrap(),
            (_, ShapeType::Solid) => model.add_solid(&kept).unwrap(),
            _ => unreachable!("a solid's own kinds"),
        };
        made.insert(shape.node(), new.clone());
        new
    }
    let bare = solid
        .located(ogeom::topo::Location::identity())
        .oriented(ogeom::topo::Orientation::Forward);
    rebuild(model, &bare, dropped.node(), &mut made)
        .located(solid.location().clone())
        .oriented(solid.orientation())
}

/// A unit cube notched along its edge at y = z = 1 by a pin of radius
/// 0.002, which takes a quarter of its disc: the notch's wall is a quarter
/// turn of a cylinder far narrower than a millimetre. Its arc at x = 0 is collapsed into one vertex whose
/// radius holds both its ends, as a file drops an edge shorter than its
/// tolerance: the wall's two rulings meet at that vertex a quarter turn
/// apart in its chart, and in space within the vertex. A box over that end
/// cuts the wall, and its boundary closes at the vertex as the topology
/// says.
#[test]
fn rulings_a_quarter_turn_apart_close_within_their_vertex() {
    let mut model = Model::with_tolerances(T);
    let radius = 0.002;
    let cube = ogeom::algo::make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
        .unwrap()
        .shape;
    let along = Frame::new(Point::new(-1.0, 1.0, 1.0), Direction::X, Direction::Y, T).unwrap();
    let pin = ogeom::algo::make_cylinder(&mut model, along, radius, 3.0, T)
        .unwrap()
        .shape;
    let notched = ogeom::boolean::cut(&mut model, &cube, &pin, T)
        .unwrap()
        .shape;
    let arc = explore_unique(&model, &notched, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let ends = explore_unique(&model, e, ShapeType::Vertex).unwrap();
            ends.len() == 2
                && ends.iter().all(|v| {
                    let p = point_of(&model, v);
                    p.x.abs() < 1e-9 && p.y > 0.99 && p.z > 0.99
                })
        })
        .expect("the notch's arc at x = 0");
    let solid = with_edge_collapsed(&mut model, &notched, &arc);
    let diagnosis = ogeom::algo::check(&model, &solid, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let notch = radius * radius * core::f64::consts::FRAC_PI_4;
    let at = Frame::new(Point::new(-0.5, 0.5, 0.5), Direction::Z, Direction::X, T).unwrap();
    let corner = ogeom::algo::make_box(&mut model, at, (1.0, 1.0, 1.0), T)
        .unwrap()
        .shape;
    let taken = 0.125 - notch / 2.0;
    for (name, want) in [("cut", 1.0 - notch - taken), ("common", taken)] {
        let made = if name == "cut" {
            ogeom::boolean::cut(&mut model, &solid, &corner, T)
        } else {
            ogeom::boolean::common(&mut model, &solid, &corner, T)
        }
        .unwrap_or_else(|e| panic!("{name}: {e}"))
        .shape;
        let diagnosis = ogeom::algo::check(&model, &made, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
        let got = volume(&model, &made);
        assert!((got - want).abs() < 1e-9, "{name}: {got} against {want}");
    }
}
