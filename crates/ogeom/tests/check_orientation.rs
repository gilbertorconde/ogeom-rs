//! What `check` sees that topology alone cannot: a face turned inside out.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use std::f64::consts::TAU;

use ogeom::algo::{
    check, make_box, make_face, make_polygon, make_revolution, make_shell, make_solid,
};
use ogeom::core::Tolerances;
use ogeom::geom::PlaneSurface;
use ogeom::math::{Axis, Direction, Frame, Plane, Point};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn cube(model: &mut Model, at: (f64, f64, f64)) -> Shape {
    let frame = Frame::new(Point::new(at.0, at.1, at.2), Direction::Z, Direction::X, T).unwrap();
    make_box(model, frame, (2.0, 3.0, 4.0), T).unwrap().shape
}

#[test]
fn a_box_is_valid_and_one_with_a_face_turned_inside_out_is_not() {
    let mut model = Model::new();
    let block = cube(&mut model, (0.0, 0.0, 0.0));
    assert!(check(&model, &block, T).unwrap().is_valid());

    let mut faces = explore_unique(&model, &block, ShapeType::Face).unwrap();
    faces[2] = faces[2].reversed();
    let shell = make_shell(&mut model, &faces).unwrap().shape;
    let flipped = make_solid(&mut model, &[shell]).unwrap().shape;
    let diagnosis = check(&model, &flipped, T).unwrap();
    assert!(!diagnosis.is_valid());
    let inward: Vec<_> = diagnosis
        .problems
        .iter()
        .filter(|p| p.kind == ShapeType::Face && p.what.contains("orientation is reversed"))
        .collect();
    assert_eq!(inward.len(), 1, "{diagnosis}");
    assert_eq!(inward[0].at.node(), faces[2].node());
    // Its ring now turns against its outward normal, so each of its four
    // edges is walked the same way by the face beside it.
    let same_way = diagnosis
        .problems
        .iter()
        .filter(|p| p.kind == ShapeType::Edge && p.what.contains("walked the same way"))
        .count();
    assert_eq!(same_way, 4, "{diagnosis}");
}

/// A drum of radius 5 and `height` along `y`, its top a cone rising
/// `rise` to the axis: a cap all but flat and not planar.
fn domed_drum(model: &mut Model, height: f64, rise: f64) -> Shape {
    let pts = [(0.0, 0.0), (5.0, 0.0), (5.0, height), (0.0, height + rise)]
        .map(|(x, y)| Point::new(x, y, 0.0));
    let wire = make_polygon(model, &pts, true, T).unwrap().shape;
    let plane = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
    let profile = make_face(model, plane, &[wire], T).unwrap().shape;
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Y,
    };
    make_revolution(model, &profile, axis, TAU, T)
        .unwrap()
        .shape
}

/// Each face of a drum with a nearly flat domed cap turned to face into
/// the drum is named, and only that face: the cap, the wall and the base,
/// on a tall drum and on one thinner than the probe's first step.
#[test]
fn a_face_of_a_drum_with_a_nearly_flat_cap_turned_inside_out_is_named() {
    for (height, rise) in [(10.0, 0.05), (10.0, 1.0), (0.05, 0.02)] {
        let mut model = Model::new();
        let drum = domed_drum(&mut model, height, rise);
        let diagnosis = check(&model, &drum, T).unwrap();
        assert!(diagnosis.is_valid(), "{height} {rise}: {diagnosis}");
        let faces = explore_unique(&model, &drum, ShapeType::Face).unwrap();
        assert_eq!(faces.len(), 3);
        for turned in 0..faces.len() {
            let mut held = faces.clone();
            held[turned] = held[turned].reversed();
            let shell = make_shell(&mut model, &held).unwrap().shape;
            let flipped = make_solid(&mut model, &[shell]).unwrap().shape;
            let diagnosis = check(&model, &flipped, T).unwrap();
            let named: Vec<_> = diagnosis
                .problems
                .iter()
                .filter(|p| p.what.contains("orientation is reversed"))
                .map(|p| p.at.node())
                .collect();
            assert_eq!(
                named,
                [faces[turned].node()],
                "{height} {rise}, face {turned}: {diagnosis}"
            );
        }
    }
}

/// A block whose top steps back over itself by under a micron: two long
/// faces a hair apart in height, joined by a sliver that faces back under
/// the lower one. Every face is turned the way its neighbours say, and the
/// block is valid. A step off the sliver longer than it is wide passes
/// under the lower face or over the upper one, and reads their sides: its
/// outside lands in the material and its inside in the air. Such a step is
/// no evidence about the sliver.
#[test]
fn a_sliver_face_narrower_than_the_probe_is_not_called_inside_out() {
    let mut model = Model::new();
    let (rise, back) = (0.0007, 0.0005);
    let pts = [
        (0.0, 0.0),
        (20.0, 0.0),
        (20.0, 10.0),
        (10.0, 10.0),
        (10.0 + back, 10.0 + rise),
        (0.0, 10.0 + rise),
    ]
    .map(|(x, y)| Point::new(x, y, 0.0));
    let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let plane = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
    let profile = make_face(&mut model, plane, &[wire], T).unwrap().shape;
    let block = ogeom::algo::make_prism(
        &mut model,
        &profile,
        ogeom::math::Vector::new(0.0, 0.0, 50.0),
        T,
    )
    .unwrap()
    .shape;
    let diagnosis = check(&model, &block, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    // The volume asks the same probe which faces to turn back, turns none,
    // and is the closed form: the profile's area times the height.
    let area = 20.0 * 10.0 + 10.0 * rise + 0.5 * back * rise;
    let volume =
        ogeom::algo::volume_properties(&model, &block, ogeom::mesh::Deflection::default(), T)
            .unwrap();
    assert_eq!(volume.deflection, 0.0, "the closed form was taken");
    let volume = volume.mass;
    assert!(
        (volume - area * 50.0).abs() < 1e-9 * area * 50.0,
        "{volume} against {}",
        area * 50.0
    );
}

/// Rigid motions to check a solid under: what `check` finds is the
/// solid's, not where it stands or which mesh triangle ranks first.
fn motions() -> Vec<ogeom::math::Transform> {
    use ogeom::math::{Transform, Vector};
    (0..24)
        .map(|step| {
            let angle = f64::from(step) * 0.262;
            Transform::translation(Vector::new(37.3, -11.9, 70.1))
                * Transform::rotation(Axis::Z, angle)
                * Transform::rotation(Axis::X, angle * 0.7)
        })
        .collect()
}

/// A prism two high over a limaçon with an inner loop, `r = 1 + 2 cos t`:
/// the profile crosses itself at the origin, so the wall swept from it
/// passes through itself along the line over the origin.
fn looped_prism(model: &mut Model) -> Shape {
    use ogeom::algo::fit::{Spacing, interpolate};
    use ogeom::geom::Curve;
    let mut pts: Vec<Point> = (0..48)
        .map(|i| {
            let t = TAU * f64::from(i) / 48.0;
            let r = 2.0f64.mul_add(t.cos(), 1.0);
            Point::new(r * t.cos(), r * t.sin(), 0.0)
        })
        .collect();
    pts.push(pts[0]);
    let spline = interpolate(&pts, 3, Spacing::Centripetal, T).unwrap();
    let domain = spline.knots().domain();
    let edge = ogeom::algo::make_edge(model, Curve::BSpline(spline), domain, T)
        .unwrap()
        .shape;
    let wire = ogeom::algo::make_wire(model, &[edge], T).unwrap().shape;
    let plane = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
    let profile = make_face(model, plane, &[wire], T).unwrap().shape;
    ogeom::algo::make_prism(model, &profile, ogeom::math::Vector::new(0.0, 0.0, 2.0), T)
        .unwrap()
        .shape
}

/// A wall that passes through itself is reported broken wherever the
/// solid stands. Probing its faces' sides cannot see it: wherever a probe
/// lands, the wall reads facing out, and the solid passed as valid.
#[test]
fn a_wall_passing_through_itself_is_broken_wherever_it_stands() {
    let mut model = Model::new();
    let prism = looped_prism(&mut model);
    let wall = explore_unique(&model, &prism, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let surface = model.node(f).unwrap().data().as_face().unwrap().surface;
            !matches!(
                model.geometry().surface(surface).unwrap(),
                ogeom::geom::SurfaceGeometry::Plane(_)
            )
        })
        .unwrap();
    for motion in motions() {
        let placed = ogeom::algo::transformed(&mut model, &prism, motion)
            .unwrap()
            .shape;
        let diagnosis = check(&model, &placed, T).unwrap();
        assert!(!diagnosis.is_usable(), "{diagnosis}");
        let crossing: Vec<_> = diagnosis
            .problems
            .iter()
            .filter(|p| p.what.contains("passes through itself"))
            .collect();
        assert_eq!(crossing.len(), 1, "{diagnosis}");
        assert_eq!(crossing[0].at.node(), wall.node(), "{diagnosis}");
    }
}

/// Two overlapping boxes as the shells of one solid. Six pairs of their
/// faces cross, each along a segment through both faces' insides (the
/// first's top and the second's side at x = 1, the first's two sides at
/// y = 0 and y = 1 with the second's bottom and its side, the first's end
/// at x = 2 with the second's bottom), and every one is named under every
/// placement. The first's top reads facing out where the second does not
/// cover it and into the material where it does, and is named for that
/// rather than settled by whichever half one probe lands on.
#[test]
fn shells_crossing_each_other_are_broken_whatever_the_placement() {
    let mut model = Model::new();
    let first = make_box(&mut model, Frame::WORLD, (2.0, 1.0, 1.0), T)
        .unwrap()
        .shape;
    let at = Frame::WORLD.with_origin(Point::new(1.0, -1.0, 0.5));
    let second = make_box(&mut model, at, (2.0, 3.0, 1.0), T).unwrap().shape;
    let shells =
        [&first, &second].map(|s| explore_unique(&model, s, ShapeType::Shell).unwrap()[0].clone());
    let solid = model.add_solid(&shells).unwrap();
    let top = explore_unique(&model, &first, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            ogeom::algo::shape_bounds(&model, f, T)
                .unwrap()
                .low()
                .unwrap()
                .z
                > 0.5
        })
        .unwrap();
    for motion in motions() {
        let placed = ogeom::algo::transformed(&mut model, &solid, motion)
            .unwrap()
            .shape;
        let diagnosis = check(&model, &placed, T).unwrap();
        let crossings = diagnosis
            .problems
            .iter()
            .filter(|p| p.what.contains("passes through another face"))
            .count();
        assert_eq!(crossings, 6, "{diagnosis}");
        assert!(
            diagnosis
                .problems
                .iter()
                .any(|p| p.at.node() == top.node() && p.what.contains("reads facing out")),
            "{diagnosis}"
        );
    }
}

/// A square turned by `turn` about the z axis at height `z`, its sides
/// `2 * half.0` and `2 * half.1`.
fn turned(model: &mut Model, z: f64, turn: f64, half: (f64, f64)) -> Shape {
    let (c, s) = (turn.cos(), turn.sin());
    let corners = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)].map(|(x, y)| {
        Point::new(
            (c * x).mul_add(half.0, -s * y * half.1),
            (s * x).mul_add(half.0, c * y * half.1),
            z,
        )
    });
    make_polygon(model, &corners, true, T).unwrap().shape
}

/// `make_solid` turns a closed shell outward by more than one probe. A
/// skin through twisted squares and rectangles carried on to a point can
/// fold near the point, its wall reading into the material there and out
/// everywhere else; one probe landing in the fold turned the whole shell
/// over. Several places agreeing settle it, and where they disagree the
/// volume the shell encloses does: given either way round, the solid made
/// encloses its volume positively.
#[test]
fn make_solid_is_not_turned_over_by_one_probe() {
    for apex_at in [
        Point::new(1.0, 2.0, 14.0),
        Point::new(0.0, 0.0, 14.0),
        Point::new(-2.0, 1.0, 12.0),
        Point::new(8.0, -6.0, 13.0),
    ] {
        let mut model = Model::new();
        let mut through = vec![
            turned(&mut model, 0.0, 0.0, (3.0, 3.0)),
            turned(&mut model, 3.0, 0.7, (6.0, 1.0)),
            turned(&mut model, 6.0, 1.5, (3.0, 3.0)),
            turned(&mut model, 9.0, 2.3, (6.0, 1.0)),
        ];
        through.push(model.add_vertex(ogeom::topo::VertexData::new(apex_at)));
        let loft = ogeom::offset::make_loft_skinned(&mut model, &through, 1e-3, T)
            .unwrap()
            .shape;
        let shell = explore_unique(&model, &loft, ShapeType::Shell).unwrap()[0].clone();
        for given in [shell.clone(), shell.reversed()] {
            let solid = make_solid(&mut model, &[given]).unwrap().shape;
            let mesh =
                ogeom::mesh::triangulate(&model, &solid, ogeom::mesh::Deflection::default(), T)
                    .unwrap();
            assert!(mesh.is_closed(), "{apex_at:?}");
            assert!(mesh.volume() > 0.0, "{apex_at:?}: {}", mesh.volume());
        }
    }
}

/// Two cubes as the shells of one solid, the second's corner pushed
/// `s = 0.01` into the first's along the diagonal: their faces cross in
/// six segments, all within `s√3` of the first's corner vertex, and are
/// named while that vertex is as tight as its edges. Once its tolerance is
/// widened past `s√3` the whole overlap lies within the vertex, which is
/// one point to its tolerance, and nothing is named.
#[test]
fn a_crossing_within_a_vertex_tolerance_is_not_named() {
    let s = 0.01;
    let mut model = Model::new();
    let first = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
        .unwrap()
        .shape;
    let at = Frame::WORLD.with_origin(Point::new(1.0 - s, 1.0 - s, 1.0 - s));
    let second = make_box(&mut model, at, (1.0, 1.0, 1.0), T).unwrap().shape;
    let shells =
        [&first, &second].map(|s| explore_unique(&model, s, ShapeType::Shell).unwrap()[0].clone());
    let solid = model.add_solid(&shells).unwrap();
    let crossings = |model: &Model| {
        check(model, &solid, T)
            .unwrap()
            .problems
            .iter()
            .filter(|p| p.what.contains("passes through another face"))
            .count()
    };
    assert_eq!(crossings(&model), 6);
    let corner = explore_unique(&model, &first, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .find(|v| {
            let at = model.node(v).unwrap().data().as_vertex().unwrap().point;
            at.distance(Point::new(1.0, 1.0, 1.0)) < 1e-9
        })
        .unwrap();
    // Twice s: past s√3, the overlap's far corner from the vertex.
    model
        .widen(&corner, ogeom::core::Tolerance::new(2.0 * s).unwrap())
        .unwrap();
    assert_eq!(crossings(&model), 0);
}

/// Sound solids whose faces meet themselves or each other without
/// crossing are valid wherever they stand: a box rounded along every
/// edge, its bands and corner patches meeting their neighbours tangent; a
/// torus, its one face meeting itself along both seams; a cone, its wall
/// closing to a point at the apex; a sphere, its face closing at both
/// poles; and a tube whose bore comes within 0.01 of its outside, where
/// the two walls' meshes cross and the surfaces do not. A crossing test
/// that took a tangent junction, a seam, a point where a surface closes,
/// or two meshes crossing for a crossing would name one of these.
#[test]
fn faces_meeting_at_a_seam_an_apex_or_tangent_are_no_crossing() {
    let mut model = Model::new();
    // The bore's centre is turned off the outside's seam by half the
    // angle the outside's mesh steps by (thirteen steps at the default
    // deflection), so the thin wall stands where the outside's chord sags
    // inward farthest, some 0.29, and the bore's mesh passes through it.
    let tube = {
        let outside = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 10.0, 5.0, T)
            .unwrap()
            .shape;
        let (off, turn) = (4.99, std::f64::consts::PI / 13.0);
        let at = Frame::WORLD.with_origin(Point::new(off * turn.cos(), off * turn.sin(), -1.0));
        let bore = ogeom::algo::make_cylinder(&mut model, at, 5.0, 7.0, T)
            .unwrap()
            .shape;
        ogeom::boolean::cut(&mut model, &outside, &bore, T)
            .unwrap()
            .shape
    };
    let block = make_box(&mut model, Frame::WORLD, (30.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let solids = [
        (
            "rounded box",
            ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
                .unwrap()
                .shape,
        ),
        (
            "torus",
            ogeom::algo::make_torus(&mut model, Frame::WORLD, 5.0, 2.0, T)
                .unwrap()
                .shape,
        ),
        (
            "cone",
            ogeom::algo::make_cone(&mut model, Frame::WORLD, 4.0, 0.0, 6.0, T)
                .unwrap()
                .shape,
        ),
        (
            "sphere",
            ogeom::algo::make_sphere(&mut model, Frame::WORLD, 3.0, T)
                .unwrap()
                .shape,
        ),
        ("thin-walled tube", tube),
    ];
    let placements = [
        ogeom::math::Transform::IDENTITY,
        motions()[5],
        motions()[17],
    ];
    for (what, solid) in solids {
        for (k, motion) in placements.iter().enumerate() {
            let placed = ogeom::algo::transformed(&mut model, &solid, *motion)
                .unwrap()
                .shape;
            let diagnosis = check(&model, &placed, T).unwrap();
            assert!(diagnosis.is_valid(), "{what}, placement {k}: {diagnosis}");
        }
    }
}
