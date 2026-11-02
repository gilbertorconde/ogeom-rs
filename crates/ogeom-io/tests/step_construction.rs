//! Construction geometry through STEP: named points, axes, trimmed curves
//! and bounded surfaces written in a geometric set beside a part's solid,
//! and read back by name, measured against the geometry written.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_doc::{Construction, ConstructionGeometry, Document};
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_geom::{BSplineCurve, CircleCurve, Curve, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Axis, Circle, Direction, Frame, KnotVector, Plane, Point, Vector};
use ogeom_topo::Model;

const T: Tolerances = Tolerances::millimetres();

/// How far a read item may stand from the one written.
const CLOSE: f64 = 1e-9;

fn tilted_frame(origin: Point) -> Frame {
    let z = Direction::new(Vector::new(0.3, -0.4, 0.866), T).unwrap();
    let x = Direction::new(Vector::new(1.0, 0.2, 0.0), T).unwrap();
    Frame::new(origin, z, x, T).unwrap()
}

/// The five items of the issue and an axis, each named.
fn items() -> Vec<(&'static str, ConstructionGeometry)> {
    let line = LineCurve::new(Axis {
        location: Point::new(1.0, 2.0, 3.0),
        direction: Direction::new(Vector::new(1.0, 1.0, 0.5), T).unwrap(),
    });
    let circle =
        CircleCurve::new(Circle::new(tilted_frame(Point::new(-4.0, 5.0, 6.0)), 7.5, T).unwrap());
    let spline = BSplineCurve::new(
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
            "edge line",
            ConstructionGeometry::Curve {
                curve: line.into(),
                range: (-2.5, 17.75),
            },
        ),
        (
            "bore rim",
            ConstructionGeometry::Curve {
                curve: circle.into(),
                range: (0.3, 2.5),
            },
        ),
        (
            "section",
            ConstructionGeometry::Curve {
                curve: spline.into(),
                range: (0.0, 1.0),
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
            "spindle axis",
            ConstructionGeometry::Axis(tilted_frame(Point::new(0.0, 0.0, 25.0))),
        ),
    ]
}

fn points_of(geometry: &ConstructionGeometry) -> Vec<Point> {
    let fractions = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0];
    match geometry {
        ConstructionGeometry::Point(p) => vec![*p],
        ConstructionGeometry::Axis(frame) => vec![
            frame.origin(),
            frame.origin() + frame.x().vector(),
            frame.origin() + frame.y().vector(),
            frame.origin() + frame.z().vector(),
        ],
        ConstructionGeometry::Curve { curve, range } => fractions
            .iter()
            .map(|f| {
                curve
                    .point_at(range.0 + f * (range.1 - range.0), T)
                    .unwrap()
            })
            .collect(),
        ConstructionGeometry::Surface { surface, u, v } => fractions
            .iter()
            .flat_map(|fu| fractions.iter().map(move |fv| (fu, fv)))
            .map(|(fu, fv)| {
                surface
                    .point_at(u.0 + fu * (u.1 - u.0), v.0 + fv * (v.1 - v.0), T)
                    .unwrap()
            })
            .collect(),
    }
}

/// The same geometry: the same kind, and the same point at every sampled
/// parameter, the stretch or window included.
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
        _ => {}
    }
    for (p, q) in points_of(written).iter().zip(points_of(read).iter()) {
        assert!(p.distance(*q) <= CLOSE, "{name}: {p:?} read as {q:?}");
    }
}

fn read_back(document: &Document) -> Document {
    let text = ogeom_io::write_step(document, T).unwrap();
    ogeom_io::read_step(&text, T).unwrap().document
}

#[test]
fn construction_beside_a_solid_reads_back_by_name() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let mut document = Document::over(model);
    let part = document.add_part("bracket", solid);
    for (name, geometry) in items() {
        document.add_construction(part, name, geometry).unwrap();
    }

    let text = ogeom_io::write_step(&document, T).unwrap();
    for keyword in [
        "GEOMETRICALLY_BOUNDED_SURFACE_SHAPE_REPRESENTATION('bracket'",
        "GEOMETRIC_SET(",
        "TRIMMED_CURVE('bore rim'",
        "RECTANGULAR_TRIMMED_SURFACE('mid plane'",
        "AXIS2_PLACEMENT_3D('spindle axis'",
    ] {
        assert!(text.contains(keyword), "{keyword} written");
    }
    let import = ogeom_io::read_step(&text, T).unwrap();
    assert_eq!(import.solids.len(), 1, "the solid still reads");
    let back = import.document;
    let (id, product) = back.products().next().unwrap();
    assert_eq!(product.name, "bracket");
    let read: Vec<&Construction> = back.construction_of(id).collect();
    assert_eq!(read.len(), items().len(), "{read:?}");
    for (name, written) in items() {
        let found = read
            .iter()
            .find(|item| item.name == name)
            .unwrap_or_else(|| panic!("{name} kept its name"));
        assert_same(name, &written, &found.geometry);
    }
}

#[test]
fn a_part_of_construction_alone_reads_back_as_one() {
    let mut model = Model::new();
    let nothing = ogeom_algo::make_compound(&mut model, &[]).unwrap().shape;
    let mut document = Document::over(model);
    let part = document.add_part("datums", nothing);
    // Points and curves alone go out as wireframe.
    let wireframe: Vec<_> = items()
        .into_iter()
        .filter(|(_, g)| !matches!(g, ConstructionGeometry::Surface { .. }))
        .collect();
    for (name, geometry) in wireframe.clone() {
        document.add_construction(part, name, geometry).unwrap();
    }

    let text = ogeom_io::write_step(&document, T).unwrap();
    assert!(text.contains("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION('datums'"));
    assert!(text.contains("GEOMETRIC_CURVE_SET("));
    let back = read_back(&document);
    let (id, product) = back.products().next().unwrap();
    assert_eq!(product.name, "datums");
    let read: Vec<&Construction> = back.construction_of(id).collect();
    assert_eq!(read.len(), wireframe.len());
    for ((name, written), found) in wireframe.iter().zip(&read) {
        assert_eq!(found.name, *name, "file order kept");
        assert_same(name, written, &found.geometry);
    }
}

#[test]
fn a_whole_circle_and_a_reversed_one_keep_their_points() {
    let circle =
        CircleCurve::new(Circle::new(tilted_frame(Point::new(1.0, 1.0, 1.0)), 3.0, T).unwrap());
    let reversed = ogeom_geom::Reversible::reversed(&Curve::from(circle));
    let mut model = Model::new();
    let nothing = ogeom_algo::make_compound(&mut model, &[]).unwrap().shape;
    let mut document = Document::over(model);
    let part = document.add_part("rims", nothing);
    let whole = ConstructionGeometry::Curve {
        curve: circle.into(),
        range: (0.0, core::f64::consts::TAU),
    };
    let turned = ConstructionGeometry::Curve {
        curve: reversed,
        range: (0.5, 2.0),
    };
    document
        .add_construction(part, "whole", whole.clone())
        .unwrap();
    document
        .add_construction(part, "turned", turned.clone())
        .unwrap();
    let back = read_back(&document);
    let (id, _) = back.products().next().unwrap();
    let read: Vec<&Construction> = back.construction_of(id).collect();
    assert_same("whole", &whole, &read[0].geometry);
    // The reversed circle's parameter is its own; its stretch is the same
    // points from the same end.
    let ConstructionGeometry::Curve { curve, range } = &read[1].geometry else {
        panic!("a curve");
    };
    for (p, f) in points_of(&turned)
        .iter()
        .zip([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])
    {
        let q = curve
            .point_at(range.0 + f * (range.1 - range.0), T)
            .unwrap();
        let on = ogeom_algo::project_on_curve(curve, *p, 64, T).unwrap();
        assert!(on.distance <= CLOSE, "turned: {p:?} is on the curve read");
        if f == 0.0 || f == 1.0 {
            assert!(p.distance(q) <= CLOSE, "turned: end {p:?} read as {q:?}");
        }
    }
}

#[test]
fn construction_persists_in_the_native_document() {
    let mut model = Model::new();
    let nothing = ogeom_algo::make_compound(&mut model, &[]).unwrap().shape;
    let mut document = Document::over(model);
    let part = document.add_part("datums", nothing);
    for (name, geometry) in items() {
        document.add_construction(part, name, geometry).unwrap();
    }
    let text =
        ogeom_io::native::write_document(&document, ogeom_io::native::WriteOptions::default())
            .unwrap();
    let back = ogeom_io::native::read_document(&text).unwrap();
    assert_eq!(back.construction(), document.construction());
}

#[test]
fn an_unknown_product_or_an_empty_range_is_refused() {
    let mut model = Model::new();
    let nothing = ogeom_algo::make_compound(&mut model, &[]).unwrap().shape;
    let mut other = Document::new();
    let foreign = other.add_part("elsewhere", nothing.clone());
    let mut document = Document::over(model);
    let line: Curve = LineCurve::new(Axis {
        location: Point::ORIGIN,
        direction: Direction::X,
    })
    .into();
    let empty = ConstructionGeometry::Curve {
        curve: line.clone(),
        range: (1.0, 1.0),
    };
    assert!(
        document
            .add_construction(foreign, "x", empty.clone())
            .is_err()
    );
    let part = document.add_part("datums", nothing);
    assert!(document.add_construction(part, "x", empty).is_err());
    let unbounded = ConstructionGeometry::Surface {
        surface: SurfaceGeometry::from(
            PlaneSurface::over(Plane::new(Frame::WORLD), (-1.0, 1.0), (-1.0, 1.0)).unwrap(),
        ),
        u: (0.0, f64::INFINITY),
        v: (0.0, 1.0),
    };
    assert!(document.add_construction(part, "y", unbounded).is_err());
}

/// A NIST test part ties a surface representation holding a geometric set
/// of two composite curves to its product; they read as the product's
/// construction curves, and the solid reads as it did.
#[test]
fn a_geometric_set_in_a_corpus_part_reads_as_construction() {
    let path = format!(
        "{}/../../tests/corpus/nist_ctc_05_asme1_rd.stp",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(path).unwrap();
    let import = ogeom_io::read_step(&text, T).unwrap();
    assert_eq!(import.solids.len(), 1);
    let curves: Vec<&Construction> = import.document.construction().iter().collect();
    assert_eq!(curves.len(), 2, "{curves:?}");
    for item in &curves {
        let ConstructionGeometry::Curve { curve, range } = &item.geometry else {
            panic!("a curve: {item:?}");
        };
        assert!(range.0 < range.1);
        assert!(curve.point_at(range.0, T).is_ok());
    }
    assert!(!import.report.skipped.contains_key("GEOMETRIC_SET"));
}
