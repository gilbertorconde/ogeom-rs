//! A face whose ring runs out along an edge and straight back: a slit, the
//! face on both sides of it, ending at a vertex no other edge reaches. It
//! bounds no material, and a boolean through the face splits the region the
//! face covers.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, edge_vertices, make_box, make_edge_between, make_face_with_pcurves, make_shell,
    make_solid, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{Curve3d as _, LineCurve, PlanarCurve, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{EdgeRepr, Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

fn block(model: &mut Model, corner: (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
    let at = Frame::new(
        Point::new(corner.0, corner.1, corner.2),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    make_box(model, at, size, T).unwrap().shape
}

/// A box 4 by 4 by 2 whose top face's ring runs from its corner at
/// (0, 0, 2) out to the middle of the top, (2, 2, 2), and straight back
/// before going round the top's four sides.
fn box_with_a_slit_in_its_top(model: &mut Model) -> Shape {
    let solid = block(model, (0.0, 0.0, 0.0), (4.0, 4.0, 2.0));
    let faces = explore_unique(model, &solid, ShapeType::Face).unwrap();
    let is_top = |model: &Model, face: &Shape| {
        let bounds = ogeom::algo::face_bounds(model, face).unwrap();
        bounds.low().unwrap().z > 1.9
    };
    let top = faces.iter().find(|f| is_top(model, f)).unwrap().clone();
    let stored = top.oriented(ogeom::topo::Orientation::Forward);
    let wire = model.children_of(&stored).unwrap()[0].clone();
    let mut ring = model.ordered_children_of(&wire).unwrap();
    let corner_at =
        |model: &Model, v: &Shape| model.node(v).unwrap().data().as_vertex().unwrap().point;
    let start = ring
        .iter()
        .position(|e| {
            let (from, _) = edge_vertices(model, e).unwrap().unwrap();
            corner_at(model, &from).distance(Point::new(0.0, 0.0, 2.0)) < 1e-9
        })
        .unwrap();
    ring.rotate_left(start);
    let (corner, _) = edge_vertices(model, &ring[0]).unwrap().unwrap();
    let (from, to) = (Point::new(0.0, 0.0, 2.0), Point::new(2.0, 2.0, 2.0));
    let tip = model.add_point(to);
    let slit = make_edge_between(
        model,
        LineCurve::segment(from, to, T).unwrap().into(),
        (0.0, from.distance(to)),
        &corner,
        &tip,
        T,
    )
    .unwrap()
    .shape;
    let mut walked = vec![slit.clone(), slit.reversed()];
    walked.extend(ring);
    let surface = model.node(&top).unwrap().data().as_face().unwrap().surface;
    let surface = model.geometry().surface(surface).unwrap().clone();
    let slitted = make_face_with_pcurves(model, surface, &[walked], T)
        .unwrap()
        .shape
        .oriented(top.orientation());
    let mut kept: Vec<Shape> = faces.into_iter().filter(|f| !f.is_same(&top)).collect();
    kept.push(slitted);
    let shell = make_shell(model, &kept).unwrap().shape;
    make_solid(model, &[shell]).unwrap().shape
}

#[test]
fn a_box_with_a_slit_in_its_top_cuts_and_meets_through_it() {
    // A block over the slit's tip and half its length, and one over the top
    // clear of the slit: each takes a cube of its own, 2 by 2 by 1 and 1 by
    // 1 by 1, out of the 32 the box holds.
    for (corner, size, taken) in [
        ((1.0, 1.0, 1.0), (2.0, 2.0, 2.0), 4.0),
        ((2.5, -1.0, 1.0), (1.0, 2.0, 2.0), 1.0),
    ] {
        let mut model = Model::with_tolerances(T);
        let solid = box_with_a_slit_in_its_top(&mut model);
        let found = check(&model, &solid, T).unwrap();
        assert!(found.is_valid(), "{found}");
        assert!((volume(&model, &solid) - 32.0).abs() <= 1e-9);
        let tool = block(&mut model, corner, size);
        for (name, want) in [("cut", 32.0 - taken), ("common", taken)] {
            let made = if name == "cut" {
                ogeom::boolean::cut(&mut model, &solid, &tool, T)
            } else {
                ogeom::boolean::common(&mut model, &solid, &tool, T)
            }
            .unwrap_or_else(|e| panic!("{name} {corner:?}: {e}"))
            .shape;
            let found = check(&model, &made, T).unwrap();
            assert!(found.is_valid(), "{name} {corner:?}: {found}");
            let got = volume(&model, &made);
            assert!(
                (got - want).abs() <= 1e-9,
                "{name} {corner:?}: {got} against {want}"
            );
        }
    }
}

/// The ring `ring` of a face with `old` replaced, where it walks it, by the
/// edges `new` (in `old`'s own direction).
fn replaced(ring: &[Shape], old: &Shape, new: &[Shape]) -> Vec<Shape> {
    let mut out = Vec::new();
    for edge in ring {
        if !edge.is_same(old) {
            out.push(edge.clone());
        } else if edge.orientation() == old.orientation() {
            out.extend(new.iter().cloned());
        } else {
            out.extend(new.iter().rev().map(Shape::reversed));
        }
    }
    out
}

/// A quarter drum of radius 2 and height 4 about the z axis, over x, y > 0,
/// whose wall's ring runs from the middle of its bottom arc, at
/// `u = pi / 4`, up the ruling there to half height and straight back. The
/// slit is stated as a seam is, its far side a period off, as a file read
/// in states an edge its face walks both ways.
fn quarter_drum_with_a_seamed_slit(model: &mut Model) -> Shape {
    let drum = ogeom::algo::make_cylinder(model, Frame::WORLD, 2.0, 4.0, T)
        .unwrap()
        .shape;
    let quadrant = block(model, (0.0, 0.0, -1.0), (3.0, 3.0, 6.0));
    let quarter = ogeom::boolean::common(model, &drum, &quadrant, T)
        .unwrap()
        .shape;
    let faces = explore_unique(model, &quarter, ShapeType::Face).unwrap();
    let surface_of = |model: &Model, face: &Shape| {
        let id = model.node(face).unwrap().data().as_face().unwrap().surface;
        model.geometry().surface(id).unwrap().clone()
    };
    let wall = faces
        .iter()
        .find(|f| matches!(surface_of(model, f), SurfaceGeometry::Cylinder(_)))
        .unwrap()
        .clone();
    let bottom = faces
        .iter()
        .find(|f| {
            ogeom::algo::face_bounds(model, f)
                .unwrap()
                .high()
                .unwrap()
                .z
                < 0.1
        })
        .unwrap()
        .clone();
    let ring_of = |model: &Model, face: &Shape| {
        let stored = face.oriented(ogeom::topo::Orientation::Forward);
        let wire = model.children_of(&stored).unwrap()[0].clone();
        model.ordered_children_of(&wire).unwrap()
    };
    let (wall_ring, bottom_ring) = (ring_of(model, &wall), ring_of(model, &bottom));
    let arc = wall_ring
        .iter()
        .find(|e| bottom_ring.iter().any(|b| b.is_same(e)))
        .unwrap()
        .clone();
    let data = model.node(&arc).unwrap().data().as_edge().unwrap().clone();
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        unreachable!("an arc has a curve")
    };
    let circle = model.geometry().curve(*curve).unwrap().clone();
    let middle = f64::midpoint(range.0, range.1);
    let at = circle.point_at(middle, T).unwrap();
    let foot = Point::new(2.0_f64.sqrt(), 2.0_f64.sqrt(), 0.0);
    assert!(at.distance(foot) < 1e-9, "{at:?}");
    let (first, last) = edge_vertices(model, &arc.oriented(ogeom::topo::Orientation::Forward))
        .unwrap()
        .unwrap();
    let split = model.add_point(foot);
    let halves = [
        make_edge_between(model, circle.clone(), (range.0, middle), &first, &split, T)
            .unwrap()
            .shape,
        make_edge_between(model, circle, (middle, range.1), &split, &last, T)
            .unwrap()
            .shape,
    ];
    let forward = arc.oriented(ogeom::topo::Orientation::Forward);
    let bottom_ring = replaced(&bottom_ring, &forward, &halves);
    let mut wall_ring = replaced(&wall_ring, &forward, &halves);
    let top = Point::new(foot.x, foot.y, 2.0);
    let tip = model.add_point(top);
    let slit = make_edge_between(
        model,
        LineCurve::segment(foot, top, T).unwrap().into(),
        (0.0, 2.0),
        &split,
        &tip,
        T,
    )
    .unwrap()
    .shape;
    // Out from the split, where the ring reaches it, and straight back.
    let reached = wall_ring
        .iter()
        .position(|e| {
            let (_, end) = edge_vertices(model, e).unwrap().unwrap();
            end.is_same(&split)
        })
        .unwrap();
    wall_ring.splice(reached + 1..reached + 1, [slit.clone(), slit.reversed()]);
    let rebuilt_wall = make_face_with_pcurves(model, surface_of(model, &wall), &[wall_ring], T)
        .unwrap()
        .shape
        .oriented(wall.orientation());
    let rebuilt_bottom =
        make_face_with_pcurves(model, surface_of(model, &bottom), &[bottom_ring], T)
            .unwrap()
            .shape
            .oriented(bottom.orientation());
    let on = model
        .node(&rebuilt_wall)
        .unwrap()
        .data()
        .as_face()
        .unwrap()
        .surface;
    let column = |model: &mut Model, u: f64| {
        let line = ogeom::geom::Line2d::segment(
            ogeom::math::Point2::new(u, 0.0),
            ogeom::math::Point2::new(u, 2.0),
            T,
        )
        .unwrap();
        model.geometry_mut().add_pcurve(PlanarCurve::Line(line))
    };
    let quarter_turn = core::f64::consts::FRAC_PI_4;
    let (near, far) = (
        column(model, quarter_turn),
        column(model, quarter_turn - core::f64::consts::TAU),
    );
    let Some(NodeData::Edge(slit_data)) = model.node_mut(&slit).map(|n| n.data_mut()) else {
        unreachable!("an edge")
    };
    slit_data
        .representations
        .retain(|r| r.surface() != Some(on));
    slit_data.representations.push(EdgeRepr::Seam {
        forward: near,
        reversed: far,
        surface: on,
        location: ogeom::topo::Location::identity(),
        range: (0.0, 2.0),
    });
    let mut kept: Vec<Shape> = faces
        .into_iter()
        .filter(|f| !f.is_same(&wall) && !f.is_same(&bottom))
        .collect();
    kept.extend([rebuilt_wall, rebuilt_bottom]);
    let shell = make_shell(model, &kept).unwrap().shape;
    make_solid(model, &[shell]).unwrap().shape
}

#[test]
fn a_quarter_drum_whose_wall_holds_a_seamed_slit_cuts_and_meets_through_it() {
    // A block over x, y > 1 and 1 < z < 3 takes the region of the quarter
    // disc past x = 1 and y = 1, of area pi / 3 - sqrt 3 + 1, over a height
    // of 2: across the slit's upper half and round its tip.
    let mut model = Model::with_tolerances(T);
    let solid = quarter_drum_with_a_seamed_slit(&mut model);
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let whole = 4.0 * core::f64::consts::PI;
    assert!((volume(&model, &solid) - whole).abs() <= 1e-6);
    let taken = 2.0 * (core::f64::consts::FRAC_PI_3 - 3.0_f64.sqrt() + 1.0);
    let tool = block(&mut model, (1.0, 1.0, 1.0), (2.0, 2.0, 2.0));
    for (name, want) in [("cut", whole - taken), ("common", taken)] {
        let made = if name == "cut" {
            ogeom::boolean::cut(&mut model, &solid, &tool, T)
        } else {
            ogeom::boolean::common(&mut model, &solid, &tool, T)
        }
        .unwrap_or_else(|e| panic!("{name}: {e}"))
        .shape;
        let found = check(&model, &made, T).unwrap();
        assert!(found.is_valid(), "{name}: {found}");
        let got = volume(&model, &made);
        assert!((got - want).abs() <= 1e-6, "{name}: {got} against {want}");
    }
}
