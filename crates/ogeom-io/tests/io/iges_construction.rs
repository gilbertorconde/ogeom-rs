//! Construction geometry through IGES: named points, axes, curves and
//! bounded surfaces written as independent entities flagged as
//! construction beside a part's solid, read back by name and measured
//! against the geometry written.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_doc::{Construction, ConstructionGeometry, Document};
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_geom::{
    BSplineCurve, BSplineSurface, CircleCurve, Curve, CylinderSurface, EllipseCurve, HelixCurve,
    HyperbolaCurve, LineCurve, ParabolaCurve, PlaneSurface, SurfaceGeometry,
};
use ogeom_math::{
    Axis, Circle, ControlGrid, Cylinder, Direction, Ellipse, Frame, Hyperbola, KnotVector,
    Parabola, Plane, Point, Vector,
};
use ogeom_topo::Model;

const T: Tolerances = Tolerances::millimetres();

/// How far a read item may stand from the one written.
const CLOSE: f64 = 1e-9;

const FRACTIONS: [f64; 7] = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0];

fn tilted_frame(origin: Point) -> Frame {
    let z = Direction::new(Vector::new(0.3, -0.4, 0.866), T).unwrap();
    let x = Direction::new(Vector::new(1.0, 0.2, 0.0), T).unwrap();
    Frame::new(origin, z, x, T).unwrap()
}

fn spline() -> BSplineCurve {
    BSplineCurve::new(
        KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0], 3).unwrap(),
        vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(10.0, 3.0, 1.0),
            Point::new(20.0, -2.0, 4.0),
            Point::new(30.0, 5.0, -1.0),
            Point::new(40.0, 0.0, 2.0),
        ],
        T,
    )
    .unwrap()
}

/// A bicubic-by-quadratic patch, not flat.
fn patch() -> BSplineSurface {
    let mut points = Vec::new();
    for i in 0..4 {
        for j in 0..3 {
            let (x, y) = (f64::from(i) * 10.0, f64::from(j) * 8.0);
            points.push(Point::new(x, y, (x * 0.3).sin() * 4.0 + y * 0.2));
        }
    }
    BSplineSurface::new(
        KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap(),
        KnotVector::new(vec![0.0, 0.0, 0.0, 2.0, 2.0, 2.0], 2).unwrap(),
        &ControlGrid::new(points, 4, 3).unwrap(),
        T,
    )
    .unwrap()
}

/// Each kind of item, each named; every one but the line reads back over
/// the parameters it was written with.
fn items() -> Vec<(&'static str, ConstructionGeometry)> {
    let circle =
        CircleCurve::new(Circle::new(tilted_frame(Point::new(-4.0, 5.0, 6.0)), 7.5, T).unwrap());
    let ellipse = EllipseCurve::new(
        Ellipse::new(tilted_frame(Point::new(3.0, -1.0, 2.0)), 9.0, 4.0, T).unwrap(),
    );
    let hyperbola = HyperbolaCurve::new(
        Hyperbola::new(tilted_frame(Point::new(1.0, 1.0, -3.0)), 3.0, 2.0, T).unwrap(),
        4.0,
    )
    .unwrap();
    let parabola = ParabolaCurve::new(
        Parabola::new(tilted_frame(Point::new(-2.0, 0.5, 1.0)), 1.25, T).unwrap(),
        20.0,
    )
    .unwrap();
    let plane = PlaneSurface::over(
        Plane::new(tilted_frame(Point::new(50.0, -5.0, 2.0))),
        (-1e3, 1e3),
        (-1e3, 1e3),
    )
    .unwrap();
    vec![
        (
            "datum point",
            ConstructionGeometry::Point(Point::new(12.5, -3.25, 0.125)),
        ),
        (
            "bore rim",
            ConstructionGeometry::Curve {
                curve: circle.into(),
                range: (0.3, 2.5),
            },
        ),
        (
            "whole rim",
            ConstructionGeometry::Curve {
                curve: circle.into(),
                range: (0.0, core::f64::consts::TAU),
            },
        ),
        (
            "cam",
            ConstructionGeometry::Curve {
                curve: ellipse.into(),
                range: (0.4, 3.0),
            },
        ),
        (
            "flare",
            ConstructionGeometry::Curve {
                curve: hyperbola.into(),
                range: (-0.5, 1.2),
            },
        ),
        (
            "dish",
            ConstructionGeometry::Curve {
                curve: parabola.into(),
                range: (-3.0, 5.0),
            },
        ),
        (
            "section",
            ConstructionGeometry::Curve {
                curve: spline().into(),
                range: (0.2, 0.8),
            },
        ),
        (
            "mid plane",
            ConstructionGeometry::Surface {
                surface: plane.into(),
                u: (-12.0, 8.0),
                v: (-3.0, 4.5),
            },
        ),
        (
            "skin",
            ConstructionGeometry::Surface {
                surface: patch().into(),
                u: (0.25, 0.75),
                v: (0.5, 2.0),
            },
        ),
        (
            "spindle axis",
            ConstructionGeometry::Axis(tilted_frame(Point::new(0.0, 0.0, 25.0))),
        ),
        (
            "mirror",
            ConstructionGeometry::Axis(tilted_frame(Point::new(4.0, -2.0, 1.0)).mirrored()),
        ),
    ]
}

fn points_of(geometry: &ConstructionGeometry) -> Vec<Point> {
    match geometry {
        ConstructionGeometry::Point(p) => vec![*p],
        ConstructionGeometry::Axis(frame) => vec![
            frame.origin(),
            frame.origin() + frame.x().vector(),
            frame.origin() + frame.y().vector(),
            frame.origin() + frame.z().vector(),
        ],
        ConstructionGeometry::Curve { curve, range } => FRACTIONS
            .iter()
            .map(|f| {
                curve
                    .point_at(range.0 + f * (range.1 - range.0), T)
                    .unwrap()
            })
            .collect(),
        ConstructionGeometry::Surface { surface, u, v } => FRACTIONS
            .iter()
            .flat_map(|fu| FRACTIONS.iter().map(move |fv| (fu, fv)))
            .map(|(fu, fv)| {
                surface
                    .point_at(u.0 + fu * (u.1 - u.0), v.0 + fv * (v.1 - v.0), T)
                    .unwrap()
            })
            .collect(),
    }
}

/// The same geometry: the same kind, the same stretch or window, and the
/// same point at every sampled parameter.
fn assert_same(name: &str, written: &ConstructionGeometry, read: &ConstructionGeometry) {
    assert_eq!(
        std::mem::discriminant(written),
        std::mem::discriminant(read),
        "{name}: {read:?}"
    );
    match (written, read) {
        (
            ConstructionGeometry::Curve { range: a, .. },
            ConstructionGeometry::Curve { range: b, .. },
        ) => assert!(
            (a.0 - b.0).abs() <= CLOSE && (a.1 - b.1).abs() <= CLOSE,
            "{name}: range {a:?} read as {b:?}"
        ),
        (
            ConstructionGeometry::Surface { u: ua, v: va, .. },
            ConstructionGeometry::Surface { u: ub, v: vb, .. },
        ) => assert!(
            [ua.0 - ub.0, ua.1 - ub.1, va.0 - vb.0, va.1 - vb.1]
                .iter()
                .all(|d| d.abs() <= CLOSE),
            "{name}: window {ua:?} {va:?} read as {ub:?} {vb:?}"
        ),
        (ConstructionGeometry::Axis(a), ConstructionGeometry::Axis(b)) => {
            assert_eq!(a.handedness(), b.handedness(), "{name}: handedness");
        }
        _ => {}
    }
    assert_same_points(name, written, read);
}

fn assert_same_points(name: &str, written: &ConstructionGeometry, read: &ConstructionGeometry) {
    for (p, q) in points_of(written).iter().zip(points_of(read).iter()) {
        assert!(p.distance(*q) <= CLOSE, "{name}: {p:?} read as {q:?}");
    }
}

/// A document whose one part holds `items` and, where `solid`, a box.
fn document_of(solid: bool, items: &[(&str, ConstructionGeometry)]) -> Document {
    let mut model = Model::new();
    let shape = if solid {
        ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape
    } else {
        ogeom_algo::make_compound(&mut model, &[]).unwrap().shape
    };
    let mut document = Document::over(model);
    let part = document.add_part("bracket", shape);
    for (name, geometry) in items {
        document
            .add_construction(part, *name, geometry.clone())
            .unwrap();
    }
    document
}

fn read_back(document: &Document) -> ogeom_io::IgesImport {
    let written = ogeom_io::write_iges_reported(document, T).unwrap();
    assert!(written.warnings.is_empty(), "{:?}", written.warnings);
    ogeom_io::read_iges(&written.text, T).unwrap()
}

#[test]
fn construction_beside_a_solid_reads_back_by_name() {
    let document = document_of(true, &items());
    let import = read_back(&document);
    assert_eq!(import.solids.len(), 1, "the solid still reads");
    assert!(import.wires.is_empty() && import.points.is_empty());
    assert!(
        import.report.warnings.is_empty(),
        "{:?}",
        import.report.warnings
    );
    assert!(
        import.report.skipped.is_empty(),
        "{:?}",
        import.report.skipped
    );
    let back = import.document;
    let (id, product) = back.products().next().unwrap();
    assert_eq!(product.name, "bracket");
    let read: Vec<&Construction> = back.construction_of(id).collect();
    assert_eq!(read.len(), items().len(), "{read:?}");
    for ((name, written), found) in items().iter().zip(&read) {
        assert_eq!(found.name, *name, "file order kept");
        assert_same(name, written, &found.geometry);
    }
}

#[test]
fn construction_alone_reads_back_as_a_part_of_it() {
    let document = document_of(false, &items());
    let import = read_back(&document);
    assert!(import.solids.is_empty() && import.sheets.is_empty());
    let back = import.document;
    let read = back.construction();
    assert_eq!(read.len(), items().len());
    for ((name, written), found) in items().iter().zip(read) {
        assert_eq!(found.name, *name);
        assert_same(name, written, &found.geometry);
    }
}

/// A line's entity states its ends and nothing of where its parameter
/// starts: it reads back as the same points from the same end, over the
/// same length.
#[test]
fn a_line_keeps_its_points_and_length() {
    let line = LineCurve::new(Axis {
        location: Point::new(1.0, 2.0, 3.0),
        direction: Direction::new(Vector::new(1.0, 1.0, 0.5), T).unwrap(),
    });
    let written = ConstructionGeometry::Curve {
        curve: line.into(),
        range: (-2.5, 17.75),
    };
    let import = read_back(&document_of(false, &[("edge line", written.clone())]));
    let found = &import.document.construction()[0];
    assert_eq!(found.name, "edge line");
    let ConstructionGeometry::Curve { range, .. } = &found.geometry else {
        panic!("a curve: {found:?}");
    };
    assert!((range.1 - range.0 - 20.25).abs() <= CLOSE, "{range:?}");
    assert_same_points("edge line", &written, &found.geometry);
}

/// A reversed curve's parameter is its own; its stretch reads back as the
/// same points from the same end.
#[test]
fn reversed_curves_keep_their_points() {
    let reversed = |curve: Curve| ogeom_geom::Reversible::reversed(&curve);
    let circle =
        CircleCurve::new(Circle::new(tilted_frame(Point::new(1.0, 1.0, 1.0)), 3.0, T).unwrap());
    let hyperbola = HyperbolaCurve::new(
        Hyperbola::new(tilted_frame(Point::ORIGIN), 2.0, 5.0, T).unwrap(),
        3.0,
    )
    .unwrap();
    let parabola = ParabolaCurve::new(
        Parabola::new(tilted_frame(Point::ORIGIN), 0.5, T).unwrap(),
        10.0,
    )
    .unwrap();
    let items = [
        (
            "turned rim",
            ConstructionGeometry::Curve {
                curve: reversed(circle.into()),
                range: (0.5, 2.0),
            },
        ),
        (
            "turned flare",
            ConstructionGeometry::Curve {
                curve: reversed(hyperbola.into()),
                range: (-1.0, 2.5),
            },
        ),
        (
            "turned dish",
            ConstructionGeometry::Curve {
                curve: reversed(parabola.into()),
                range: (-4.0, 6.0),
            },
        ),
    ];
    let import = read_back(&document_of(false, &items));
    let read = import.document.construction();
    assert_eq!(read.len(), items.len());
    for ((name, written), found) in items.iter().zip(read) {
        assert_eq!(found.name, *name);
        assert_same_points(name, written, &found.geometry);
    }
}

/// A cylinder's window goes out as its exact spline patch, whose parameter
/// is not the angle: it reads back as that patch over its own domain, the
/// window's corners where they were.
#[test]
fn a_curved_window_reads_back_as_its_exact_patch() {
    let cylinder = CylinderSurface::new(
        Cylinder::new(tilted_frame(Point::new(5.0, 5.0, 0.0)), 4.0, T).unwrap(),
        (-50.0, 50.0),
    )
    .unwrap();
    let (u, v) = ((0.2, 2.0), (-3.0, 7.0));
    let written = ConstructionGeometry::Surface {
        surface: cylinder.into(),
        u,
        v,
    };
    let import = read_back(&document_of(false, &[("wall", written)]));
    let found = &import.document.construction()[0];
    assert_eq!(found.name, "wall");
    let ConstructionGeometry::Surface {
        surface,
        u: ru,
        v: rv,
    } = &found.geometry
    else {
        panic!("a surface: {found:?}");
    };
    let window = ogeom_geom::TrimmedSurface::new(cylinder.into(), u, v, T).unwrap();
    let exact = SurfaceGeometry::Trimmed(Box::new(window))
        .to_bspline(T)
        .unwrap();
    assert_same_points(
        "wall",
        &ConstructionGeometry::Surface {
            surface: exact.clone().into(),
            u: exact.domain().0,
            v: exact.domain().1,
        },
        &found.geometry,
    );
    let original = SurfaceGeometry::from(cylinder);
    for (s, t, rs, rt) in [
        (u.0, v.0, ru.0, rv.0),
        (u.1, v.0, ru.1, rv.0),
        (u.0, v.1, ru.0, rv.1),
        (u.1, v.1, ru.1, rv.1),
    ] {
        let p = original.point_at(s, t, T).unwrap();
        let q = surface.point_at(rs, rt, T).unwrap();
        assert!(p.distance(q) <= CLOSE, "corner {p:?} read as {q:?}");
    }
}

/// A helix has no exact spline form: it is left out with a warning naming
/// it, and the rest is written.
#[test]
fn an_item_with_no_iges_form_is_left_out_with_a_warning() {
    let helix = HelixCurve::new(Frame::WORLD, 5.0, 2.0, 3.0).unwrap();
    let document = document_of(
        true,
        &[
            (
                "datum",
                ConstructionGeometry::Point(Point::new(1.0, 2.0, 3.0)),
            ),
            (
                "thread",
                ConstructionGeometry::Curve {
                    curve: helix.into(),
                    range: (0.0, 6.0),
                },
            ),
            ("origin", ConstructionGeometry::Axis(Frame::WORLD)),
        ],
    );
    let written = ogeom_io::write_iges_reported(&document, T).unwrap();
    assert_eq!(written.warnings.len(), 1, "{:?}", written.warnings);
    assert!(
        written.warnings[0].contains("'thread'"),
        "{:?}",
        written.warnings
    );
    let import = ogeom_io::read_iges(&written.text, T).unwrap();
    assert!(
        import.report.skipped.is_empty(),
        "{:?}",
        import.report.skipped
    );
    let names: Vec<&str> = import
        .document
        .construction()
        .iter()
        .map(|c| c.name.as_str())
        .collect();
    assert_eq!(names, ["datum", "origin"]);
}

/// A name a Name property cannot carry whole goes out as its label, and
/// the writer says so.
#[test]
fn a_name_with_a_delimiter_keeps_its_label_with_a_warning() {
    let document = document_of(
        false,
        &[("pin, left", ConstructionGeometry::Point(Point::ORIGIN))],
    );
    let written = ogeom_io::write_iges_reported(&document, T).unwrap();
    assert_eq!(written.warnings.len(), 1, "{:?}", written.warnings);
    assert!(
        written.warnings[0].contains("'pin, lef'"),
        "{:?}",
        written.warnings
    );
    let import = ogeom_io::read_iges(&written.text, T).unwrap();
    assert_eq!(import.document.construction()[0].name, "pin, lef");
}

/// A document with no construction geometry writes as it did: no entity
/// carries the construction flag and nothing is reported.
#[test]
fn a_document_without_construction_writes_no_construction() {
    let document = document_of(true, &[]);
    let written = ogeom_io::write_iges_reported(&document, T).unwrap();
    assert!(written.warnings.is_empty());
    assert_eq!(written.text, ogeom_io::write_iges(&document, T).unwrap());
    for record in written
        .text
        .lines()
        .filter(|l| l.as_bytes().get(72) == Some(&b'D'))
    {
        let status = &record[64..72];
        let odd_line = record[73..].trim().parse::<u32>().unwrap() % 2 == 1;
        assert!(!odd_line || &status[4..6] == "00", "{record}");
    }
    let import = ogeom_io::read_iges(&written.text, T).unwrap();
    assert!(import.document.construction().is_empty());
}
