//! A thread's groove: an M6 tooth space swept seven turns down a helix,
//! and cut from a block.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::Containment;
use ogeom::algo::{
    check, classify_in_solid_exact, make_box, make_cylinder, make_edge, make_face, make_polygon,
    make_wire, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{Curve2d as _, Curve3d as _, HelixCurve, PlaneSurface, Surface as _};
use ogeom::math::{Direction, Frame, Plane, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{EdgeRepr, Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

/// The groove: radial from 2.4 to 3.0 about an axis at (10, 10) running
/// down from z = 11, 0.8 wide at 2.4 and 0.125 at 3.0, pitch 1, its
/// section square to the helix where it starts.
fn groove(model: &mut Model, turns: f64) -> Shape {
    let (pitch, radius) = (1.0, 2.7);
    let frame = Frame::new(Point::new(10.0, 10.0, 11.0), -Direction::Z, Direction::X, T).unwrap();
    let helix = HelixCurve::new(frame, radius, pitch, turns).unwrap();
    let range = helix.domain();
    let edge = make_edge(model, helix.into(), range, T).unwrap().shape;
    let spine = make_wire(model, &[edge], T).unwrap().shape;
    let x = frame.x().vector();
    let tangent = (frame.y().vector() * radius
        + frame.z().vector() * (pitch / core::f64::consts::TAU))
        .normalized(T)
        .unwrap();
    let across = tangent.cross(x).normalized(T).unwrap();
    let axis = frame.origin();
    let corners: Vec<Point> = [(2.4, -0.4), (3.0, -0.0625), (3.0, 0.0625), (2.4, 0.4)]
        .iter()
        .map(|&(r, w)| axis + x * r + across * w)
        .collect();
    let wire = make_polygon(model, &corners, true, T).unwrap().shape;
    let plane = Plane::through(axis + x * radius, Direction::new(tangent, T).unwrap());
    let section = make_face(model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape;
    ogeom::offset::make_pipe_shell(model, &section, &spine, false, 1e-3, T)
        .unwrap()
        .shape
}

/// Every edge of the groove lies where each face it bounds says it does:
/// its image on the face, read through the face's surface, is on the edge's
/// curve at the same parameter, to the edge's tolerance. Two walls sharing a
/// rail fit their surfaces at their own paces along the helix, and an image
/// assumed straight on the wall that adopted the rail drifted a fifth of a
/// millimetre off it.
#[test]
fn a_swept_groove_is_where_its_faces_say() {
    let mut model = Model::new();
    let groove = groove(&mut model, 7.0);
    assert!(check(&model, &groove, T).unwrap().is_valid());
    for face in explore_unique(&model, &groove, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        let surface = model.geometry().surface(data.surface).unwrap().clone();
        for edge in explore_unique(&model, &face, ShapeType::Edge).unwrap() {
            let ed = model.node(&edge).unwrap().data().as_edge().unwrap().clone();
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = ed.curve3d() else {
                continue;
            };
            let curve = model.geometry().curve(*curve).unwrap().clone();
            let Some(EdgeRepr::PCurve {
                curve: image,
                range: image_range,
                ..
            }) = ed.pcurve_for(data.surface, edge.location())
            else {
                continue;
            };
            let image = model.geometry().pcurve(*image).unwrap().clone();
            for k in 0..=32 {
                let f = f64::from(k) / 32.0;
                let t = range.0 + (range.1 - range.0) * f;
                let s = image_range.0 + (image_range.1 - image_range.0) * f;
                let uv = image.point_at(s, T).unwrap();
                let ((u0, u1), (v0, v1)) = surface.domain();
                let off = surface
                    .point_at(uv.x.clamp(u0, u1), uv.y.clamp(v0, v1), T)
                    .unwrap()
                    .distance(curve.point_at(t, T).unwrap());
                assert!(
                    off <= ed.tolerance.get() + 1e-6,
                    "an edge's image stands {off:e} off it, against {:e}",
                    ed.tolerance.get()
                );
            }
        }
    }
}

/// The groove cut from a block for one, two and seven turns: valid, and
/// the block less exactly what the groove holds of it.
#[test]
fn a_thread_groove_cuts_from_a_block() {
    for turns in [1.0, 2.0, 7.0] {
        let mut model = Model::new();
        let groove = groove(&mut model, turns);
        let block = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &block, &groove, T)
            .unwrap_or_else(|e| panic!("{turns} turns: {e}"))
            .shape;
        let diagnosis = check(&model, &cut, T).unwrap();
        assert!(diagnosis.is_valid(), "{turns} turns: {diagnosis}");
        let inside = ogeom::boolean::common(&mut model, &block, &groove, T)
            .unwrap()
            .shape;
        let (left, taken) = (volume(&model, &cut), volume(&model, &inside));
        assert!(
            (left + taken - 4000.0).abs() < 4000.0 * 1e-6,
            "{turns} turns: {left} + {taken}"
        );
        assert!(taken > 3.0 * (turns - 1.0), "{turns} turns take {taken}");
    }
}

/// A point a hundredth of a millimetre behind the end of a three-turn
/// groove, near its wall: inside. A ray from it leaves through the wall a
/// few hundredths on, close to where the wall's patch ends, and the Newton
/// solve seeded at the corner of the patch's last cell stepped past the
/// patch's end and stalled against it, so the crossing went uncounted.
#[test]
fn a_point_just_inside_the_grooves_end_is_inside() {
    let mut model = Model::new();
    let groove = groove(&mut model, 3.0);
    let point = Point::new(
        12.499_990_690_040_166,
        9.993_177_265_028_299,
        8.231_476_863_123_2,
    );
    assert_eq!(
        classify_in_solid_exact(&model, &groove, point, T).unwrap(),
        Containment::In
    );
}

/// The groove cut from a block bored 2.5 in radius up the groove's axis
/// from z = 2, so the groove's inner flank runs in and out of the bore's
/// wall: valid, and the bored block less what the groove holds of it.
fn cuts_from_a_bored_block(turns: f64) {
    let mut model = Model::new();
    let groove = groove(&mut model, turns);
    let block = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let axis = Frame::new(Point::new(10.0, 10.0, 2.0), Direction::Z, Direction::X, T).unwrap();
    let bore = make_cylinder(&mut model, axis, 2.5, 9.0, T).unwrap().shape;
    let bored = ogeom::boolean::cut(&mut model, &block, &bore, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut model, &bored, &groove, T)
        .unwrap_or_else(|e| panic!("{turns} turns: {e}"))
        .shape;
    let diagnosis = check(&model, &cut, T).unwrap();
    assert!(diagnosis.is_valid(), "{turns} turns: {diagnosis}");
    let inside = ogeom::boolean::common(&mut model, &bored, &groove, T)
        .unwrap_or_else(|e| panic!("{turns} turns: {e}"))
        .shape;
    let diagnosis = check(&model, &inside, T).unwrap();
    assert!(diagnosis.is_valid(), "{turns} turns: {diagnosis}");
    // The bore's wall, cut into pieces, is measured by its mesh, and the
    // mesh of a hollow lies inside it by a share of the chord over its
    // whole area. At a chord of 3e-5 the two add up to a millionth.
    let chord = 1e-3;
    let whole = 4000.0 - core::f64::consts::PI * 2.5 * 2.5 * 8.0;
    let wall = core::f64::consts::TAU * 2.5 * 8.0;
    let (left, taken) = (volume(&model, &cut), volume(&model, &inside));
    assert!(
        (left + taken - whole).abs() < whole * 1e-6 + wall * chord,
        "{turns} turns: {left} + {taken} against {whole}"
    );
    assert!(taken > 3.0 * (turns - 1.0), "{turns} turns take {taken}");
}

/// Three turns end the groove's cap on the line where the bore's surface
/// closes on itself, leaving a sliver of the bore's wall beside that line
/// inside the groove.
#[test]
fn a_three_turn_groove_cuts_from_a_bored_block() {
    cuts_from_a_bored_block(3.0);
}

/// Six turns end it there too, where the groove's rail grazes the bore:
/// the sections on the two walls either side of the rail meet it a
/// hundredth of a millimetre apart.
#[test]
fn a_six_turn_groove_cuts_from_a_bored_block() {
    cuts_from_a_bored_block(6.0);
}

/// Seven turns end the cap across the bore's wall, and the flank's section
/// winds round the bore for turns before it.
#[test]
fn a_seven_turn_groove_cuts_from_a_bored_block() {
    cuts_from_a_bored_block(7.0);
}

/// Five and three quarter turns end the cap where one flank, turned about
/// the helix by the sweep, has just dipped inside the bore: the flank meets
/// the bore's wall along a two-millimetre run at a grazing angle, too
/// shallow for sampled cells to cross, found from where it leaves the
/// flank's border.
#[test]
fn a_groove_whose_flank_just_dips_into_the_bore_cuts_from_it() {
    cuts_from_a_bored_block(5.75);
    cuts_from_a_bored_block(6.75);
}
