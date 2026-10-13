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
