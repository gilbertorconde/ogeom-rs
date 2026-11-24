//! The `.brep` interchange format: what the writer emits, the reader
//! returns, and a file written by hand from the specification reads as the
//! shape it describes.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

fn counts(model: &Model, shape: &Shape) -> (usize, usize, usize, usize) {
    (
        explore_unique(model, shape, ShapeType::Vertex)
            .unwrap()
            .len(),
        explore_unique(model, shape, ShapeType::Edge).unwrap().len(),
        explore_unique(model, shape, ShapeType::Face).unwrap().len(),
        explore_unique(model, shape, ShapeType::Solid)
            .unwrap()
            .len(),
    )
}

#[test]
fn a_drilled_block_survives_the_round_trip() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 12.0), T)
        .unwrap()
        .shape;
    let bore = ogeom::algo::make_cylinder(
        &mut model,
        Frame::new(Point::new(30.0, 15.0, -1.0), Direction::Z, Direction::X, T).unwrap(),
        4.0,
        14.0,
        T,
    )
    .unwrap()
    .shape;
    let part = ogeom::boolean::cut(&mut model, &block, &bore, T)
        .unwrap()
        .shape;

    let text = ogeom::io::brep::write(&model, &part, T).unwrap();
    assert!(
        text.starts_with("DBRep_DrawableShape"),
        "the file says what it is"
    );

    let (back, shape) = ogeom::io::brep::read(&text, T).unwrap();
    assert_eq!(
        counts(&back, &shape),
        counts(&model, &part),
        "every vertex, edge, face and solid comes back"
    );
    let (was, now) = (volume(&model, &part), volume(&back, &shape));
    assert!(
        (now - was).abs() < was * 1e-9,
        "and it is the same solid: {now} against {was}"
    );

    // And the text is a fixed point from there: write, read, write again,
    // and the bytes are identical, so the format carries everything the
    // reader put back and the reader takes everything the format carries.
    //
    // The comparison starts at the *second* text rather than the first,
    // because one flag legitimately changes on the way in: the boolean's
    // own edges never had their same-parameter claim re-established, and
    // the reader measures the representations and establishes it.
    let again = ogeom::io::brep::write(&back, &shape, T).unwrap();
    let (twice, shape) = ogeom::io::brep::read(&again, T).unwrap();
    let third = ogeom::io::brep::write(&twice, &shape, T).unwrap();
    assert_eq!(third, again, "the text is a fixed point");
}

#[test]
fn seams_and_poles_come_back_as_themselves() {
    // A drum and a ball: between them a cylinder with a seam, a sphere with
    // a seam and two degenerate poles, two planes, and circular edges. The
    // format has to carry all of it, including the edges that are points.
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T)
        .unwrap()
        .shape;
    let ball = ogeom::algo::make_sphere(
        &mut model,
        Frame::new(Point::new(30.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
        5.0,
        T,
    )
    .unwrap()
    .shape;
    let both = ogeom::algo::make_compound(&mut model, &[drum, ball])
        .unwrap()
        .shape;

    let text = ogeom::io::brep::write(&model, &both, T).unwrap();
    let (back, shape) = ogeom::io::brep::read(&text, T).unwrap();
    assert_eq!(counts(&back, &shape), counts(&model, &both));

    // The sphere is still a sphere and the cylinder still a cylinder, read
    // back as the surfaces they are, not fitted into splines.
    let kinds = |model: &Model, shape: &Shape| -> (usize, usize) {
        let faces = explore_unique(model, shape, ShapeType::Face).unwrap();
        let mut spheres = 0;
        let mut cylinders = 0;
        for face in faces {
            let ogeom::topo::NodeData::Face(data) = model.node(&face).unwrap().data() else {
                continue;
            };
            match model.geometry().surface(data.surface) {
                Some(ogeom::geom::SurfaceGeometry::Sphere(_)) => spheres += 1,
                Some(ogeom::geom::SurfaceGeometry::Cylinder(_)) => cylinders += 1,
                _ => {}
            }
        }
        (spheres, cylinders)
    };
    assert_eq!(kinds(&back, &shape), kinds(&model, &both));
    assert_eq!(kinds(&back, &shape), (1, 1));

    // And the degenerate edges (the sphere's poles) survived as degenerate
    // edges rather than as nothing.
    let poles = |model: &Model, shape: &Shape| -> usize {
        explore_unique(model, shape, ShapeType::Edge)
            .unwrap()
            .iter()
            .filter(|e| {
                model
                    .node(e)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| d.degenerate)
            })
            .count()
    };
    assert_eq!(poles(&back, &shape), poles(&model, &both));
    assert_eq!(poles(&back, &shape), 2, "the ball has two poles");
}

#[test]
fn a_file_written_by_hand_reads_as_the_square_it_describes() {
    // Four vertices, four edges on four lines, one wire, one planar face,
    // spelled out against the specification rather than produced by the
    // writer, so the reader is tested against the format and not against
    // its own dialect.
    let text = "\
DBRep_DrawableShape

CASCADE Topology V1, (c) Matra-Datavision
Locations 0
Curve2ds 4
1 0 0 1 0
1 10 0 0 1
1 10 10 -1 0
1 0 10 0 -1
Curves 4
1 0 0 0 1 0 0
1 10 0 0 0 1 0
1 10 10 0 -1 0 0
1 0 10 0 0 -1 0
Polygon3D 0
PolygonOnTriangulations 0
Surfaces 1
1 0 0 0 0 0 1 1 0 0 0 1 0
Triangulations 0

TShapes 10
Ve
1e-07
0 0 0
0 0

0001000
*
Ve
1e-07
10 0 0
0 0

0001000
*
Ed
 1e-07 1 1 0
1 1 0 0 10
2 1 1 0 0 10
0

0001000
+10 0 -9 0 *
Ve
1e-07
10 10 0
0 0

0001000
*
Ed
 1e-07 1 1 0
1 2 0 0 10
2 2 1 0 0 10
0

0001000
+9 0 -7 0 *
Ve
1e-07
0 10 0
0 0

0001000
*
Ed
 1e-07 1 1 0
1 3 0 0 10
2 3 1 0 0 10
0

0001000
+7 0 -5 0 *
Ed
 1e-07 1 1 0
1 4 0 0 10
2 4 1 0 0 10
0

0001000
+5 0 -10 0 *
Wi

0001000
+8 0 +6 0 +4 0 +3 0 *
Fa
0 1e-07 1 0

0001000
+2 0 *

+1 0
";
    let (model, face) = ogeom::io::brep::read(text, T).unwrap();
    assert_eq!(model.kind_of(&face).unwrap(), ShapeType::Face);
    let (vertices, edges, faces, _) = counts(&model, &face);
    assert_eq!((vertices, edges, faces), (4, 4, 1));

    // It is the ten-by-ten square in the z = 0 plane: its area says so.
    let area = ogeom::algo::surface_properties(&model, &face, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((area - 100.0).abs() < 1e-9, "a ten by ten square: {area}");
}

#[test]
fn placed_shapes_and_cone_charts_survive_the_round_trip() {
    // A cone and an L-shaped prism, the pair placed by a turn and a shift:
    // the placement has to be written as a placement record the reader
    // takes back, the cone's charts in the format's own terms (along the
    // ruling), and the prism's edges charted where each face uses them.
    let mut model = Model::new();
    let cone = ogeom::algo::make_cone(&mut model, Frame::WORLD, 5.0, 2.0, 10.0, T)
        .unwrap()
        .shape;
    let corners = [
        Point::new(20.0, 0.0, 0.0),
        Point::new(40.0, 0.0, 0.0),
        Point::new(40.0, 5.0, 0.0),
        Point::new(25.0, 5.0, 0.0),
        Point::new(25.0, 15.0, 0.0),
        Point::new(20.0, 15.0, 0.0),
    ];
    let outline = ogeom::algo::make_polygon(&mut model, &corners, true, T)
        .unwrap()
        .shape;
    let plane = ogeom::geom::PlaneSurface::new(ogeom::math::Plane::new(Frame::WORLD)).into();
    let ell = ogeom::algo::make_face(&mut model, plane, &[outline], T)
        .unwrap()
        .shape;
    let prism =
        ogeom::algo::make_prism(&mut model, &ell, ogeom::math::Vector::new(0.0, 0.0, 7.0), T)
            .unwrap()
            .shape;
    let both = ogeom::algo::make_compound(&mut model, &[cone, prism])
        .unwrap()
        .shape;
    let turn = ogeom::math::Transform::rotation(
        ogeom::math::Axis::new(Point::new(1.0, 2.0, 3.0), Direction::X),
        0.7,
    ) * ogeom::math::Transform::translation(ogeom::math::Vector::new(5.0, -4.0, 2.0));
    let datum = model.add_datum(turn);
    let placed = both.moved(&ogeom::topo::Location::of(datum));

    let text = ogeom::io::brep::write(&model, &placed, T).unwrap();
    let (back, shape) = ogeom::io::brep::read(&text, T).unwrap();
    assert_eq!(counts(&back, &shape), counts(&model, &placed));
    // A frustum of radii 5 and 2 over 10, and an L of area 150 over 7: both
    // integrate in closed form, so a billionth is rounding's share.
    let frustum = core::f64::consts::PI * 10.0 / 3.0 * (25.0 + 10.0 + 4.0);
    let exact = frustum + 150.0 * 7.0;
    let now = ogeom::algo::volume_properties(&back, &shape, Deflection::default(), T).unwrap();
    assert!(
        (now.mass - exact).abs() < exact * 1e-9,
        "{} against {exact}",
        now.mass
    );
    let was = ogeom::algo::volume_properties(&model, &placed, Deflection::default(), T).unwrap();
    assert!(
        now.centre.distance(was.centre) < 1e-9,
        "and where it was: {:?} against {:?}",
        now.centre,
        was.centre
    );
}

#[test]
fn a_read_extrusion_survives_the_native_format() {
    // An elliptic prism's wall is an extrusion; read back from `.brep` it
    // reaches both ways from its curve, and the native format has to carry
    // that window through a write and a read.
    let mut model = Model::new();
    let ellipse = ogeom::geom::Curve::Ellipse(ogeom::geom::EllipseCurve::new(
        ogeom::math::Ellipse::new(Frame::WORLD, 10.0, 5.0, T).unwrap(),
    ));
    let range = ogeom::geom::Curve3d::domain(&ellipse);
    let edge = ogeom::algo::make_edge(&mut model, ellipse, range, T)
        .unwrap()
        .shape;
    let plane =
        ogeom::geom::PlaneSurface::over(ogeom::math::Plane::XY, (-20.0, 20.0), (-20.0, 20.0))
            .unwrap();
    let face = ogeom::algo::make_face_with_pcurves(&mut model, plane.into(), &[vec![edge]], T)
        .unwrap()
        .shape;
    let prism = ogeom::algo::make_prism(
        &mut model,
        &face,
        ogeom::math::Vector::new(0.0, 0.0, 4.0),
        T,
    )
    .unwrap()
    .shape;

    let text = ogeom::io::brep::write(&model, &prism, T).unwrap();
    let (read, solid) = ogeom::io::brep::read(&text, T).unwrap();
    assert!(
        read.geometry()
            .surfaces()
            .any(|(_, s)| matches!(s, ogeom::geom::SurfaceGeometry::Extrusion(_))),
        "the wall is an extrusion"
    );
    let options = ogeom::io::native::WriteOptions::default();
    let native = ogeom::io::native::write(&read, std::slice::from_ref(&solid), options).unwrap();
    let (back, roots) = ogeom::io::native::read(&native).unwrap();
    assert_eq!(
        ogeom::io::native::write(&back, &roots, options).unwrap(),
        native
    );
    assert_eq!(counts(&back, &roots[0]), counts(&read, &solid));
    assert!(ogeom::algo::check(&back, &roots[0], T).unwrap().is_valid());
    let (now, was) = (volume(&back, &roots[0]), volume(&model, &prism));
    assert!((now - was).abs() < 1e-9 * was, "{now} against {was}");
}

/// A curve in the `z = 0` plane paired with its chart curve on that plane
/// under the world chart, at the same parameter, and the part an edge uses.
type Drawn = (
    &'static str,
    ogeom::geom::Curve,
    ogeom::geom::PlanarCurve,
    (f64, f64),
);

/// One of each curve kind the format has a record for, forward, each with
/// its chart curve on the plane: lines, conics, splines, and trims and
/// offsets of them, a trim across a circle's seam and one of a curve
/// already run backwards among them.
fn drawn_curves() -> Vec<Drawn> {
    use ogeom::geom::{
        BSpline2d, BSplineCurve, Circle2d, CircleCurve, Curve, Curve3d as _, Ellipse2d,
        EllipseCurve, Hyperbola2d, HyperbolaCurve, Line2d, LineCurve, Offset2d, OffsetCurve,
        ParabolaCurve, PlanarCurve, Reversible as _, Trimmed2d, TrimmedCurve,
    };
    use ogeom::math::{
        Axis, Axis2, Circle, Circle2, Direction2, Ellipse, Ellipse2, Frame2, Hyperbola, Hyperbola2,
        KnotVector, Parabola, Point2, Vector, Vector2,
    };
    let (angle, centre) = (0.4_f64, Point2::new(0.3, -0.2));
    let x = Vector::new(angle.cos(), angle.sin(), 0.0);
    let frame = Frame::new(
        Point::new(centre.x, centre.y, 0.0),
        Direction::Z,
        Direction::new(x, T).unwrap(),
        T,
    )
    .unwrap();
    let frame2 = Frame2::new(centre, Direction2::new(Vector2::new(x.x, x.y), T).unwrap());
    let flat = |p: Point| Point2::new(p.x, p.y);

    let line: Curve = LineCurve::over(Axis::new(frame.origin(), frame.x()), -1.0, 2.0)
        .unwrap()
        .into();
    let line2: PlanarCurve = Line2d::over(Axis2::new(centre, frame2.x()), -1.0, 2.0)
        .unwrap()
        .into();
    let circle: Curve = CircleCurve::new(Circle::new(frame, 1.5, T).unwrap()).into();
    let circle2: PlanarCurve = Circle2d::new(Circle2::new(frame2, 1.5, T).unwrap()).into();
    let ellipse: Curve = EllipseCurve::new(Ellipse::new(frame, 2.0, 1.0, T).unwrap()).into();
    let ellipse2: PlanarCurve = Ellipse2d::new(Ellipse2::new(frame2, 2.0, 1.0, T).unwrap()).into();
    let hyperbola: Curve =
        HyperbolaCurve::over(Hyperbola::new(frame, 1.2, 0.7, T).unwrap(), -1.0, 2.0)
            .unwrap()
            .into();
    let hyperbola2: PlanarCurve =
        Hyperbola2d::over(Hyperbola2::new(frame2, 1.2, 0.7, T).unwrap(), -1.0, 2.0)
            .unwrap()
            .into();
    // A parabola is quadratic in its parameter: in the plane, the Bezier
    // through its ends whose middle pole stands on both end tangents.
    let parabola: Curve = ParabolaCurve::over(Parabola::new(frame, 0.8, T).unwrap(), -1.0, 2.0)
        .unwrap()
        .into();
    let (start, end) = (
        parabola.point_at(-1.0, T).unwrap(),
        parabola.point_at(2.0, T).unwrap(),
    );
    let middle = start + parabola.d1_at(-1.0, T).unwrap() * 1.5;
    let parabola2: PlanarCurve = BSpline2d::new(
        KnotVector::new(vec![-1.0, -1.0, -1.0, 2.0, 2.0, 2.0], 2).unwrap(),
        vec![flat(start), flat(middle), flat(end)],
        T,
    )
    .unwrap()
    .into();
    let poles = [(0.0, 0.0), (0.5, 1.0), (1.2, -0.4), (2.0, 0.8), (2.6, 0.1)];
    let knots = || KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0], 3).unwrap();
    let spline: Curve = BSplineCurve::new(
        knots(),
        poles.iter().map(|&(x, y)| Point::new(x, y, 0.0)).collect(),
        T,
    )
    .unwrap()
    .into();
    let spline2: PlanarCurve = BSpline2d::new(
        knots(),
        poles.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
        T,
    )
    .unwrap()
    .into();

    let trim = |c: &Curve, s: f64, e: f64| -> Curve {
        TrimmedCurve::new(c.clone(), s, e, T).unwrap().into()
    };
    let trim2 = |c: &PlanarCurve, s: f64, e: f64| -> PlanarCurve {
        Trimmed2d::new(c.clone(), s, e, T).unwrap().into()
    };
    // An offset along the tangent crossed with `-z` is one along the
    // chart curve's left normal.
    let offset = |c: &Curve| -> Curve {
        Curve::Offset(Box::new(
            OffsetCurve::new(c.clone(), 0.25, -Direction::Z).unwrap(),
        ))
    };
    let offset2 = |c: &PlanarCurve| -> PlanarCurve {
        PlanarCurve::Offset(Box::new(Offset2d::new(c.clone(), 0.25).unwrap()))
    };
    vec![
        ("line", line.clone(), line2.clone(), (-0.5, 1.5)),
        ("circle", circle.clone(), circle2.clone(), (0.4, 2.5)),
        ("ellipse", ellipse.clone(), ellipse2.clone(), (0.4, 2.5)),
        (
            "hyperbola",
            hyperbola.clone(),
            hyperbola2.clone(),
            (-0.5, 1.5),
        ),
        ("parabola", parabola, parabola2, (-0.5, 1.5)),
        ("spline", spline.clone(), spline2.clone(), (0.1, 0.9)),
        (
            "trimmed line",
            trim(&line, -0.8, 1.7),
            trim2(&line2, -0.8, 1.7),
            (-0.5, 1.5),
        ),
        (
            "circle trimmed across its seam",
            trim(&circle, 5.0, 7.0),
            trim2(&circle2, 5.0, 7.0),
            (5.2, 6.8),
        ),
        (
            "trimmed ellipse",
            trim(&ellipse, 0.5, 2.5),
            trim2(&ellipse2, 0.5, 2.5),
            (0.7, 2.3),
        ),
        (
            "trimmed hyperbola",
            trim(&hyperbola, -0.5, 1.8),
            trim2(&hyperbola2, -0.5, 1.8),
            (-0.2, 1.5),
        ),
        (
            "trimmed spline",
            trim(&spline, 0.2, 0.7),
            trim2(&spline2, 0.2, 0.7),
            (0.3, 0.6),
        ),
        (
            "trim of a circle run backwards",
            trim(&circle.reversed(), 1.0, 3.0),
            trim2(&circle2.reversed(), 1.0, 3.0),
            (1.2, 2.8),
        ),
        (
            "offset circle",
            offset(&circle),
            offset2(&circle2),
            (0.4, 2.5),
        ),
        (
            "offset of a trimmed hyperbola",
            offset(&trim(&hyperbola, -0.5, 1.8)),
            offset2(&trim2(&hyperbola2, -0.5, 1.8)),
            (-0.2, 1.5),
        ),
    ]
}

/// An edge on `curve` over `range` between `ends`, carrying `chart` on
/// `surface` over the same range.
fn edge_with_chart(
    model: &mut Model,
    (curve, chart): (ogeom::geom::Curve, ogeom::geom::PlanarCurve),
    range: (f64, f64),
    surface: ogeom::topo::SurfaceId,
    ends: &[Shape],
) -> Shape {
    use ogeom::topo::{EdgeData, EdgeRepr, Location};
    let curve = model.geometry_mut().add_curve(curve);
    let chart = model.geometry_mut().add_pcurve(chart);
    let mut data = EdgeData::on_curve(curve, Location::identity(), range);
    data.add(EdgeRepr::PCurve {
        curve: chart,
        surface,
        location: Location::identity(),
        range,
    });
    data.assert_same_parameter(true);
    model.add_edge(data, ends).unwrap()
}

/// Every curve kind crosses as an edge's curve in space and as its chart
/// curve on a plane, forward and run backwards. The edge bounds a face on
/// the plane, closed through a point off its chord. Read back, the curve
/// and the chart curve each stand where they were written at every
/// fraction of the edge's range, to 1e-12, and still share one range, so
/// they still agree at each parameter.
#[test]
fn every_curve_kind_crosses_forward_and_run_backwards() {
    use ogeom::geom::{
        Curve2d as _, Curve3d as _, Line2d, LineCurve, Reversible as _, Surface as _,
    };
    use ogeom::math::{Point2, Vector};
    use ogeom::topo::{EdgeRepr, FaceData, Location};
    let plane: ogeom::geom::SurfaceGeometry =
        ogeom::geom::PlaneSurface::new(ogeom::math::Plane::new(Frame::WORLD)).into();
    let flat = |p: Point| Point2::new(p.x, p.y);
    let mut failures = Vec::new();
    for (name, curve, chart, range) in drawn_curves() {
        assert_eq!(curve.domain(), chart.domain(), "{name}");
        let (lo, hi) = curve.domain();
        let backwards = (lo + hi - range.1, lo + hi - range.0);
        for (way, curve, chart, range) in [
            ("forward", curve.clone(), chart.clone(), range),
            ("backwards", curve.reversed(), chart.reversed(), backwards),
        ] {
            let at = format!("{name}, {way}");
            let mut model = Model::new();
            let surface = model.geometry_mut().add_surface(plane.clone());
            let (start, end) = (
                curve.point_at(range.0, T).unwrap(),
                curve.point_at(range.1, T).unwrap(),
            );
            let aside = start.midpoint(end) + Vector::new(0.0, 0.0, 1.0).cross(end - start);
            let corners = [start, end, aside].map(|p| model.add_point(p));
            let mut edges = vec![edge_with_chart(
                &mut model,
                (curve.clone(), chart.clone()),
                range,
                surface,
                &corners[..2],
            )];
            for (from, to) in [(1, 2), (2, 0)] {
                let (a, b) = ([start, end, aside][from], [start, end, aside][to]);
                let segment = (
                    LineCurve::segment(a, b, T).unwrap().into(),
                    Line2d::segment(flat(a), flat(b), T).unwrap().into(),
                );
                let edge = edge_with_chart(
                    &mut model,
                    segment,
                    (0.0, a.distance(b)),
                    surface,
                    &[corners[from].clone(), corners[to].clone()],
                );
                edges.push(edge);
            }
            let wire = model.add_wire(&edges).unwrap();
            let face = model
                .add_face(FaceData::new(surface, Location::identity()), &[wire])
                .unwrap();

            let text = match ogeom::io::brep::write(&model, &face, T) {
                Ok(text) => text,
                Err(e) => {
                    failures.push(format!("{at}: refused: {e}"));
                    continue;
                }
            };
            let (back, shape) = ogeom::io::brep::read(&text, T).unwrap();
            // The edge is the one whose ends are the curve's, read within a
            // hundredth of the confusion distance: each closing segment has
            // an end at the point aside, more than a chord off both.
            let Some(data) = explore_unique(&back, &shape, ShapeType::Edge)
                .unwrap()
                .into_iter()
                .filter_map(|e| back.node(&e).and_then(|n| n.data().as_edge()).cloned())
                .find(|data| {
                    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                        return false;
                    };
                    let read = back.geometry().curve(*curve).unwrap();
                    [(range.0, start), (range.1, end)]
                        .iter()
                        .all(|(t, p)| read.point_at(*t, T).unwrap().distance(*p) < 1e-9)
                })
            else {
                failures.push(format!("{at}: the edge is lost"));
                continue;
            };
            let Some(EdgeRepr::Curve3d {
                curve: read,
                range: read_range,
                ..
            }) = data.curve3d()
            else {
                unreachable!("found by its curve");
            };
            let Some((read_chart, chart_range, on)) =
                data.representations.iter().find_map(|r| match r {
                    EdgeRepr::PCurve {
                        curve,
                        range,
                        surface,
                        ..
                    } => Some((*curve, *range, *surface)),
                    _ => None,
                })
            else {
                failures.push(format!("{at}: the chart curve is lost"));
                continue;
            };
            let read = back.geometry().curve(*read).unwrap();
            let read_chart = back.geometry().pcurve(read_chart).unwrap();
            let on = back.geometry().surface(on).unwrap();
            if *read_range != chart_range {
                failures.push(format!(
                    "{at}: curve over {read_range:?}, chart curve over {chart_range:?}"
                ));
            }
            let mut widest = (0.0_f64, 0.0_f64);
            for k in 0..=64 {
                let f = f64::from(k) / 64.0;
                let along = |(a, b): (f64, f64)| (b - a).mul_add(f, a);
                let point = curve.point_at(along(range), T).unwrap();
                let drawn = chart.point_at(along(range), T).unwrap();
                let read_point = read.point_at(along(*read_range), T).unwrap();
                let read_drawn = read_chart.point_at(along(chart_range), T).unwrap();
                let lifted = on.point_at(read_drawn.x, read_drawn.y, T).unwrap();
                widest.0 = widest.0.max(read_point.distance(point));
                widest.1 = widest
                    .1
                    .max(read_drawn.distance(drawn))
                    .max(lifted.distance(point));
            }
            // Each curve written is the one in the model with its frame
            // turned, its origin or knots moved or its window mirrored, all
            // by arithmetic exact to rounding, on curves a few units across.
            if widest.0 > 1e-12 || widest.1 > 1e-12 {
                failures.push(format!(
                    "{at}: the curve misses by {:.3e}, the chart curve by {:.3e}",
                    widest.0, widest.1
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
