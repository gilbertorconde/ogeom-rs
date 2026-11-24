//! The middle path of pipe-like solids, pinned against spines known in
//! closed form.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::Curve3d as _;
use ogeom_math::{Circle, Frame, Point};
use ogeom_topo::{EdgeRepr, Filter, Model, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

/// The path's tolerance: a hundredth of a millimetre on parts a few
/// millimetres across.
const TOLERANCE: f64 = 0.01;

fn face_centres(model: &Model, solid: &Shape) -> Vec<(Shape, Point)> {
    let deflection = ogeom_mesh::Deflection {
        chord: 1e-3,
        ..ogeom_mesh::Deflection::default()
    };
    explore(model, solid, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .map(|f| {
            let centre = ogeom_algo::surface_properties(model, &f, deflection, T)
                .unwrap()
                .centre;
            (f, centre)
        })
        .collect()
}

/// The face whose centroid is nearest `at`.
fn face_at(model: &Model, solid: &Shape, at: Point) -> Shape {
    face_centres(model, solid)
        .into_iter()
        .min_by(|a, b| a.1.distance(at).total_cmp(&b.1.distance(at)))
        .unwrap()
        .0
}

/// The path's single edge, as a curve and its range.
fn path_curve(model: &Model, wire: &Shape) -> (ogeom_geom::Curve, (f64, f64)) {
    let edges = explore(model, wire, Filter::OfType(ShapeType::Edge)).unwrap();
    assert_eq!(edges.len(), 1);
    let data = model.node(&edges[0]).unwrap().data().as_edge().unwrap();
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        panic!("the path edge carries a 3D curve");
    };
    (model.geometry().curve(*curve).unwrap().clone(), *range)
}

/// The largest distance from points along the path to `distance_to`.
fn worst(model: &Model, wire: &Shape, distance_to: impl Fn(Point) -> f64) -> (f64, Point, Point) {
    let (curve, (lo, hi)) = path_curve(model, wire);
    let mut worst = 0.0_f64;
    for i in 0..=200 {
        let p = curve
            .point_at(lo + (hi - lo) * f64::from(i) / 200.0, T)
            .unwrap();
        worst = worst.max(distance_to(p));
    }
    (
        worst,
        curve.point_at(lo, T).unwrap(),
        curve.point_at(hi, T).unwrap(),
    )
}

#[test]
fn a_cylinders_middle_path_is_its_axis() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 1.5, 6.0, T)
        .unwrap()
        .shape;
    let bottom = face_at(&model, &solid, Point::ORIGIN);
    let top = face_at(&model, &solid, Point::new(0.0, 0.0, 6.0));

    let path = ogeom_offset::middle_path(&mut model, &solid, &bottom, &top, TOLERANCE, T).unwrap();
    assert!(path.deviation <= TOLERANCE);
    let (curve, _) = path_curve(&model, &path.built.shape);
    assert!(
        matches!(curve, ogeom_geom::Curve::Line(_)),
        "a straight tube's path is a line"
    );
    let (off, a, b) = worst(&model, &path.built.shape, |p| p.x.hypot(p.y));
    assert!(off < TOLERANCE, "{off} off the axis");
    assert!(a.distance(Point::ORIGIN) < TOLERANCE, "starts at {a:?}");
    assert!(
        b.distance(Point::new(0.0, 0.0, 6.0)) < TOLERANCE,
        "ends at {b:?}"
    );
    assert!(!path.built.history.generated(&bottom).is_empty());
    assert!(!path.built.history.generated(&top).is_empty());
}

#[test]
fn a_bored_cylinders_middle_path_is_its_axis() {
    let mut model = Model::new();
    let outer = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T)
        .unwrap()
        .shape;
    // The bore sits off the axis, so the section's centroid does too: the
    // hole must be subtracted, not ignored, for the path to land where the
    // measured centroid is.
    let bore_frame = Frame::new(
        Point::new(0.8, 0.0, -1.0),
        ogeom_math::Direction::Z,
        ogeom_math::Direction::X,
        T,
    )
    .unwrap();
    let bore = ogeom_algo::make_cylinder(&mut model, bore_frame, 0.6, 7.0, T)
        .unwrap()
        .shape;
    let solid = ogeom_bool::cut(&mut model, &outer, &bore, T).unwrap().shape;
    let solid = explore(&model, &solid, Filter::OfType(ShapeType::Solid))
        .unwrap()
        .remove(0);
    let bottom = face_at(&model, &solid, Point::ORIGIN);
    let top = face_at(&model, &solid, Point::new(0.0, 0.0, 5.0));

    let path = ogeom_offset::middle_path(&mut model, &solid, &bottom, &top, TOLERANCE, T).unwrap();
    // The annulus's centroid: the disc's at the origin less the bore's.
    let (a_out, a_in) = (4.0, 0.36);
    let x = -(0.8 * a_in) / (a_out - a_in);
    let (off, a, _) = worst(&model, &path.built.shape, |p| (p.x - x).hypot(p.y));
    assert!(
        off < TOLERANCE,
        "{off} off the section centroids' line x = {x}"
    );
    assert!((a.z).abs() < TOLERANCE);
}

#[test]
fn a_quarter_arc_pipes_middle_path_is_its_arc() {
    let mut model = Model::new();
    let circle = Circle::new(Frame::WORLD, 2.0, T).unwrap();
    let curve = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(circle));
    let spine = ogeom_algo::make_edge(&mut model, curve, (0.0, core::f64::consts::FRAC_PI_2), T)
        .unwrap()
        .shape;
    let solid = ogeom_offset::make_pipe(&mut model, &spine, 0.3, T)
        .unwrap()
        .shape;
    let start = face_at(&model, &solid, Point::new(2.0, 0.0, 0.0));
    let end = face_at(&model, &solid, Point::new(0.0, 2.0, 0.0));

    let path = ogeom_offset::middle_path(&mut model, &solid, &start, &end, TOLERANCE, T).unwrap();
    let (off, a, b) = worst(&model, &path.built.shape, |p| {
        (p.x.hypot(p.y) - 2.0).hypot(p.z)
    });
    assert!(off < TOLERANCE, "{off} off the spine arc");
    assert!(a.distance(Point::new(2.0, 0.0, 0.0)) < TOLERANCE);
    assert!(b.distance(Point::new(0.0, 2.0, 0.0)) < TOLERANCE);
}

#[test]
fn a_ring_bent_almost_shut_is_walked_the_long_way() {
    let mut model = Model::new();
    let circle = Circle::new(Frame::WORLD, 2.0, T).unwrap();
    let curve = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(circle));
    let sweep = 350.0_f64.to_radians();
    let spine = ogeom_algo::make_edge(&mut model, curve, (0.0, sweep), T)
        .unwrap()
        .shape;
    let solid = ogeom_offset::make_pipe(&mut model, &spine, 0.3, T)
        .unwrap()
        .shape;
    // The end faces stand a third of a millimetre apart, closer than the
    // march's first step.
    let (first, last) = (
        Point::new(2.0, 0.0, 0.0),
        Point::new(2.0 * sweep.cos(), 2.0 * sweep.sin(), 0.0),
    );
    let start = face_at(&model, &solid, first);
    let end = face_at(&model, &solid, last);

    let path = ogeom_offset::middle_path(&mut model, &solid, &start, &end, TOLERANCE, T).unwrap();
    let (off, a, b) = worst(&model, &path.built.shape, |p| {
        (p.x.hypot(p.y) - 2.0).hypot(p.z)
    });
    assert!(off < TOLERANCE, "{off} off the spine arc");
    assert!(a.distance(first) < TOLERANCE);
    assert!(b.distance(last) < TOLERANCE);
}

#[test]
fn a_helical_pipes_middle_path_is_its_helix() {
    let mut model = Model::new();
    let helix = ogeom_geom::HelixCurve::new(Frame::WORLD, 3.0, 4.0, 1.25).unwrap();
    let curve: ogeom_geom::Curve = helix.into();
    let domain = curve.domain();
    let samples: Vec<Point> = (0..=4000)
        .map(|i| {
            let t = domain.0 + (domain.1 - domain.0) * f64::from(i) / 4000.0;
            curve.point_at(t, T).unwrap()
        })
        .collect();
    let (first, last) = (samples[0], samples[samples.len() - 1]);
    let spine = ogeom_algo::make_edge(&mut model, curve, domain, T)
        .unwrap()
        .shape;
    // The skin to a tenth of the path's tolerance, so the tube's own fit
    // cannot account for what the path is held to.
    let solid = ogeom_offset::make_pipe_skinned(&mut model, &spine, 0.4, 1e-3, T)
        .unwrap()
        .shape;
    let start = face_at(&model, &solid, first);
    let end = face_at(&model, &solid, last);

    let path = ogeom_offset::middle_path(&mut model, &solid, &start, &end, TOLERANCE, T).unwrap();
    // Distance to the helix through its dense polyline, whose chords stand
    // off the helix by well under a micron at this spacing.
    let (off, a, b) = worst(&model, &path.built.shape, |p| {
        samples
            .windows(2)
            .map(|w| {
                let d = w[1] - w[0];
                let s = ((p - w[0]).dot(d) / d.dot(d)).clamp(0.0, 1.0);
                p.distance(w[0] + d * s)
            })
            .fold(f64::INFINITY, f64::min)
    });
    assert!(off < TOLERANCE, "{off} off the helix");
    assert!(a.distance(first) < TOLERANCE);
    assert!(b.distance(last) < TOLERANCE);
}

#[test]
fn a_middle_path_needs_two_faces_of_the_solid() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 2.0, T)
        .unwrap()
        .shape;
    let other = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 2.0, T)
        .unwrap()
        .shape;
    let bottom = face_at(&model, &solid, Point::ORIGIN);
    let foreign = face_at(&model, &other, Point::new(0.0, 0.0, 2.0));
    assert!(ogeom_offset::middle_path(&mut model, &solid, &bottom, &bottom, TOLERANCE, T).is_err());
    assert!(
        ogeom_offset::middle_path(&mut model, &solid, &bottom, &foreign, TOLERANCE, T).is_err()
    );
    assert!(ogeom_offset::middle_path(&mut model, &solid, &bottom, &foreign, -1.0, T).is_err());
}

/// A round rod of radius 2 swept along the polyline through `corners`,
/// each corner mitred.
fn mitred_rod(model: &mut Model, corners: &[Point]) -> Shape {
    use ogeom_geom::{CircleCurve, PlaneSurface};
    use ogeom_math::{Direction, Plane};
    let frame = Frame::about(
        corners[0],
        Direction::new(corners[1] - corners[0], T).unwrap(),
    );
    let circle = Circle::new(frame, 2.0, T).unwrap();
    let edge = ogeom_algo::make_edge(
        model,
        CircleCurve::new(circle).into(),
        (0.0, 2.0 * core::f64::consts::PI),
        T,
    )
    .unwrap()
    .shape;
    let wire = ogeom_algo::make_wire(model, &[edge], T).unwrap().shape;
    let disc = ogeom_algo::make_face(
        model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let spine = ogeom_algo::make_polygon(model, corners, false, T)
        .unwrap()
        .shape;
    let rod = ogeom_offset::make_pipe_shell_with(
        model,
        &disc,
        &spine,
        &ogeom_offset::PipeLaw::RotationMinimizing,
        ogeom_offset::PipeCorners::Mitre,
        1e-3,
        T,
    )
    .unwrap()
    .shape;
    assert!(ogeom_algo::check(model, &rod, T).unwrap().is_valid());
    rod
}

/// The distance from `p` to the polyline through `corners`.
fn off_polyline(corners: &[Point], p: Point) -> f64 {
    corners
        .windows(2)
        .map(|w| {
            let d = w[1] - w[0];
            let s = ((p - w[0]).dot(d) / d.dot(d)).clamp(0.0, 1.0);
            p.distance(w[0] + d * s)
        })
        .fold(f64::INFINITY, f64::min)
}

/// The middle path of a rod with one mitred corner is the rod's spine: a
/// straight edge along each leg, meeting at the corner, every point of it
/// on the polyline and its length the polyline's.
fn the_middle_path_follows_the_mitred_spine(corners: [Point; 3]) {
    let mut model = Model::new();
    let solid = mitred_rod(&mut model, &corners);
    let start = face_at(&model, &solid, corners[0]);
    let end = face_at(&model, &solid, corners[2]);
    let tolerance = 0.02;
    let path = ogeom_offset::middle_path(&mut model, &solid, &start, &end, tolerance, T).unwrap();
    assert!(path.deviation <= tolerance);

    // The wire is open, which the check takes for a gap in a boundary, so
    // its edges are checked one by one.
    let edges = explore(&model, &path.built.shape, Filter::OfType(ShapeType::Edge)).unwrap();
    assert_eq!(edges.len(), 2);
    for edge in &edges {
        let diagnosis = ogeom_algo::check(&model, edge, T).unwrap();
        assert!(diagnosis.is_valid(), "{diagnosis}");
    }
    let mut length = 0.0;
    let mut ends = Vec::new();
    let mut worst = 0.0_f64;
    for edge in &edges {
        let data = model.node(edge).unwrap().data().as_edge().unwrap();
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            panic!("a path edge carries a 3D curve");
        };
        let curve = model.geometry().curve(*curve).unwrap();
        let (lo, hi) = *range;
        let mut previous = curve.point_at(lo, T).unwrap();
        ends.push(previous);
        for i in 1..=100 {
            let p = curve
                .point_at(lo + (hi - lo) * f64::from(i) / 100.0, T)
                .unwrap();
            worst = worst.max(off_polyline(&corners, p));
            length += p.distance(previous);
            previous = p;
        }
        ends.push(previous);
    }
    assert!(worst < tolerance, "{worst} off the spine");
    let want = corners[0].distance(corners[1]) + corners[1].distance(corners[2]);
    assert!(
        (length - want).abs() < 2.0 * tolerance,
        "a length of {length} against {want}"
    );
    assert!(ends[0].distance(corners[0]) < tolerance);
    assert!(
        ends[1].distance(corners[1]) < tolerance,
        "a corner at {:?}",
        ends[1]
    );
    assert!(ends[2].distance(corners[1]) < tolerance);
    assert!(ends[3].distance(corners[2]) < tolerance);
}

#[test]
fn a_rod_with_a_mitred_square_corner_has_its_spine_for_a_middle_path() {
    the_middle_path_follows_the_mitred_spine([
        Point::ORIGIN,
        Point::new(0.0, 0.0, 20.0),
        Point::new(15.0, 0.0, 20.0),
    ]);
}

#[test]
fn rods_turning_gently_and_sharply_have_their_spines_for_middle_paths() {
    for turn in [50.0_f64, 120.0] {
        let (s, c) = turn.to_radians().sin_cos();
        the_middle_path_follows_the_mitred_spine([
            Point::ORIGIN,
            Point::new(0.0, 0.0, 20.0),
            Point::new(15.0 * s, 0.0, 20.0 + 15.0 * c),
        ]);
    }
}

/// Two sharp corners leave a leg neither march reaches: the legs from the
/// two ends do not meet, and the path is refused rather than guessed.
#[test]
fn a_rod_with_two_mitred_corners_is_refused() {
    let corners = [
        Point::ORIGIN,
        Point::new(0.0, 0.0, 20.0),
        Point::new(15.0, 0.0, 20.0),
        Point::new(15.0, 0.0, 40.0),
    ];
    let mut model = Model::new();
    let solid = mitred_rod(&mut model, &corners);
    let start = face_at(&model, &solid, corners[0]);
    let end = face_at(&model, &solid, corners[3]);
    let refused = ogeom_offset::middle_path(&mut model, &solid, &start, &end, 0.02, T).unwrap_err();
    assert!(
        refused.to_string().contains("turned back on itself"),
        "{refused}"
    );
}
