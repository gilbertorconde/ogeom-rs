//! Features whose ends lie on the body's own faces: pads stacked end to end
//! and refined into one, and a loft between sections drawn facing each
//! other, cut through a body whose end faces its sections lie on.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{Containment, make_box, make_edge, make_face, make_polygon, make_prism};
use ogeom::core::Tolerances;
use ogeom::geom::CircleCurve;
use ogeom::math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

fn assert_valid(model: &Model, shape: &Shape, what: &str) {
    let diagnosis = ogeom::algo::check(model, shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{what}: {diagnosis}");
}

/// A stepped outline drawn on the plane `y = at` facing down the y axis,
/// clockwise as seen from that side, padded `depth` further down it: each
/// wall's plane is windowed about the wall alone, and the solid holds its
/// walls reversed.
fn stepped_pad(model: &mut Model, at: f64, depth: f64) -> Shape {
    let outline = [
        (20.0, 10.0),
        (20.0, -10.0),
        (15.0, -10.0),
        (15.0, -18.0),
        (-15.0, -18.0),
        (-15.0, -10.0),
        (-20.0, -10.0),
        (-20.0, 10.0),
    ]
    .map(|(x, z)| Point::new(x, at, z));
    let wire = make_polygon(model, &outline, true, T).unwrap().shape;
    let down = Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap();
    let frame = Frame::new(Point::new(0.0, at, 0.0), down, Direction::X, T).unwrap();
    let plane = ogeom::geom::PlaneSurface::new(Plane::new(frame));
    let face = make_face(model, plane.into(), &[wire], T).unwrap().shape;
    make_prism(model, &face, Vector::new(0.0, -depth, 0.0), T)
        .unwrap()
        .shape
}

/// Two pads end to end, refined: each wall of the first merges with the
/// second's below it. The merged walls keep their material inside, walked
/// the way the solid walks them, and a point inside the second pad, beside
/// a wall whose plane was windowed about the first, still reads inside.
#[test]
fn walls_refined_across_stacked_pads_keep_the_material_inside() {
    let mut model = Model::new();
    let first = stepped_pad(&mut model, 0.0, 4.0);
    let second = stepped_pad(&mut model, -4.0, 20.0);
    let fused = ogeom::boolean::fuse(&mut model, &first, &second, T)
        .unwrap()
        .shape;
    let (refined, merged) = ogeom::heal::unify_same_domain(&mut model, &fused, T).unwrap();
    assert!(merged > 0, "the walls merge");
    assert_valid(&model, &refined.shape, "the refined pads");
    // The outline's area, 40 by 20 and 30 by 8, times the 24 mm of both
    // pads. Plane faces measure in closed form.
    let exact = (40.0 * 20.0 + 30.0 * 8.0) * 24.0;
    assert!((volume(&model, &refined.shape) - exact).abs() <= 1e-6 * exact);
    for p in [
        Point::new(19.5, -20.0, 0.0),
        Point::new(-14.5, -20.0, -17.5),
        Point::new(0.0, -23.5, 9.5),
    ] {
        assert_eq!(
            ogeom::algo::classify_in_solid_exact(&model, &refined.shape, p, T).unwrap(),
            Containment::In,
            "{p:?}"
        );
    }
}

/// A circle of radius 6.1 on the plane `y = at` facing `facing`, its centre
/// `rise` up the z axis and its start on the plane's x axis `start`.
fn section(model: &mut Model, at: f64, rise: f64, facing: f64, start: f64) -> Shape {
    let normal = Direction::new(Vector::new(0.0, facing, 0.0), T).unwrap();
    let x = Direction::new(Vector::new(start, 0.0, 0.0), T).unwrap();
    let frame = Frame::new(Point::new(0.0, at, rise), normal, x, T).unwrap();
    let curve = CircleCurve::new(Circle::new(frame, 6.1, T).unwrap());
    let edge = make_edge(model, curve.into(), (0.0, 2.0 * PI), T)
        .unwrap()
        .shape;
    ogeom::algo::make_wire(model, &[edge], T).unwrap().shape
}

/// A slanted bore lofted between circles sketched on the two end faces of
/// a block, each sketch facing out of its own face, so the two circles
/// wind opposite ways round. The loft is the slanted cylinder through
/// them, not a skin joining each point to one across the loop, and the
/// cut through the block, its ends on the block's end faces, leaves the
/// block less the cylinder.
#[test]
fn a_loft_between_sections_facing_each_other_bores_through_the_block() {
    let mut model = Model::new();
    let block = make_box(
        &mut model,
        Frame::WORLD.with_origin(Point::new(-10.0, -5.0, -10.0)),
        (20.0, 25.0, 20.0),
        T,
    )
    .unwrap()
    .shape;
    let top = section(&mut model, 20.0, -1.05, 1.0, -1.0);
    let bottom = section(&mut model, -5.0, 1.05, -1.0, 1.0);
    let tool = ogeom::offset::make_loft_skinned(&mut model, &[top, bottom], 1e-3, T)
        .unwrap()
        .shape;
    assert_valid(&model, &tool, "the loft");
    // A slanted cylinder holds its section times its height between the
    // section planes, 25 mm. The skin is fitted within 1e-3 mm, which
    // moves the volume by at most that times its 38 mm round and 25 mm
    // long wall.
    let cylinder = PI * 6.1 * 6.1 * 25.0;
    let tool_volume = volume(&model, &tool);
    assert!(
        (tool_volume - cylinder).abs() <= 1e-3 * 38.4 * 25.0,
        "{tool_volume} against {cylinder}"
    );

    let cut = ogeom::boolean::cut(&mut model, &block, &tool, T)
        .unwrap()
        .shape;
    let common = ogeom::boolean::common(&mut model, &block, &tool, T)
        .unwrap()
        .shape;
    assert_valid(&model, &cut, "the cut");
    assert_valid(&model, &common, "the common");
    // The bore lies within the block, its ends on the block's ends: what
    // is kept in common is the loft, and the two parts make the block. The
    // parts share the loft's faces, so they measure alike to rounding.
    let (cut_volume, common_volume) = (volume(&model, &cut), volume(&model, &common));
    assert!((common_volume - tool_volume).abs() <= 1e-6 * tool_volume);
    assert!((cut_volume + common_volume - 10_000.0).abs() <= 1e-6 * 10_000.0);
}
