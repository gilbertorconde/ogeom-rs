//! IGES, held to the same standard as every other format here: what goes out
//! comes back, and the measurement says so.
//!
//! Each case builds a solid with the kernel's own operations, writes it as an
//! IGES manifold solid B-rep, reads the deck back, and measures the recovered
//! solid's volume against the original's: not against a hope, against the
//! number. The corpus covers the surface vocabulary a real file exercises:
//! planes, a periodic cylinder wall with its seam, a torus, a boolean result
//! carrying both, and a B-spline loft.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::SurfaceKind;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// Write one solid, read it back, and return (original volume, recovered
/// volume, the import) for the assertions each case cares about.
fn round_trip(mut model: Model, solid: Shape) -> (f64, f64, ogeom::io::IgesImport) {
    let original = volume(&model, &solid);
    let mut document = ogeom::doc::Document::over(std::mem::take(&mut model));
    document.add_part("part", solid);
    let text = ogeom::io::write_iges(&document, T).unwrap();
    let import = ogeom::io::read_iges(&text, T).unwrap();
    assert_eq!(import.solids.len(), 1, "one solid out, one back");
    let recovered = volume(import.document.model(), &import.solids[0]);
    (original, recovered, import)
}

#[test]
fn a_box_round_trips_to_its_own_volume() {
    let mut model = Model::new();
    let solid = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 12.0, 8.0), T)
        .unwrap()
        .shape;
    let (original, recovered, import) = round_trip(model, solid);
    assert!((original - 20.0 * 12.0 * 8.0).abs() < 1e-9);
    assert!(
        (recovered - original).abs() < 1e-9,
        "{recovered} against {original}"
    );
    assert!(
        import.report.warnings.is_empty(),
        "{:?}",
        import.report.warnings
    );
}

#[test]
fn a_cylinder_round_trips_with_its_seam() {
    let mut model = Model::new();
    let solid = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 6.0, 15.0, T)
        .unwrap()
        .shape;
    let (original, recovered, _) = round_trip(model, solid);
    // The wall wraps its whole period. A reader that lost the seam would
    // refuse to tessellate at all, and one that miscounted the turn would
    // miss by a factor, not an epsilon.
    assert!(
        (recovered - original).abs() < original * 1e-6,
        "{recovered} against {original}"
    );
}

#[test]
fn a_torus_round_trips_doubly_periodic() {
    let mut model = Model::new();
    let solid = ogeom::algo::make_torus(&mut model, Frame::WORLD, 10.0, 3.0, T)
        .unwrap()
        .shape;
    let (original, recovered, _) = round_trip(model, solid);
    assert!(
        (recovered - original).abs() < original * 1e-6,
        "{recovered} against {original}"
    );
}

#[test]
fn a_drilled_block_round_trips_through_its_boolean_faces() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (30.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(15.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, frame, 4.0, 12.0, T)
        .unwrap()
        .shape;
    let solid = ogeom::boolean::cut(&mut model, &block, &drill, T)
        .unwrap()
        .shape;
    let (original, recovered, _) = round_trip(model, solid);
    // The wall's rims come back chained (each
    // image shifted by whole periods until its start meets the last one's
    // end), so the chart is the rectangle it was, and the two solids
    // tessellate the same way. A tenth of a percent between them at this
    // chord is the mesh reading a torn chart.
    assert!(
        (recovered - original).abs() < original * 1e-9,
        "{recovered} against {original}"
    );
}

#[test]
fn a_spline_walled_prism_round_trips_through_its_tabulated_cylinders() {
    // A lens-shaped profile from two fitted open B-spline arcs, extruded:
    // the walls are extrusion surfaces over B-spline curves, which the
    // writer spells as tabulated cylinders (122) over rational B-spline
    // curves (126), read back as the same extrusions, and the fitted-pcurve
    // path carries back in.
    //
    // The walls are open here. The closed-wall case has its own test below,
    // for both formats at once.
    let mut model = Model::new();
    let arc = |sign: f64| -> ogeom::geom::Curve {
        let pts: Vec<Point> = (0..=8)
            .map(|i| {
                let t = f64::from(i) / 8.0;
                let x = -10.0 + 20.0 * t;
                let y = sign * 6.0 * (core::f64::consts::PI * t).sin();
                Point::new(x, y, 0.0)
            })
            .collect();
        ogeom::geom::Curve::from(
            ogeom::geom::fit::fit_points(&pts, 3, 1e-9, T)
                .unwrap()
                .curve,
        )
    };
    let upper = arc(1.0);
    let lower = arc(-1.0);
    let va = ogeom::algo::make_vertex(&mut model, Point::new(-10.0, 0.0, 0.0)).shape;
    let vb = ogeom::algo::make_vertex(&mut model, Point::new(10.0, 0.0, 0.0)).shape;
    let (u0, u1) = ogeom::geom::Curve3d::domain(&upper);
    let (l0, l1) = ogeom::geom::Curve3d::domain(&lower);
    let e_up = ogeom::algo::make_edge_between(&mut model, upper, (u0, u1), &va, &vb, T)
        .unwrap()
        .shape;
    let e_dn = ogeom::algo::make_edge_between(&mut model, lower, (l0, l1), &va, &vb, T)
        .unwrap()
        .shape;
    let face = ogeom::algo::make_face_with_pcurves(
        &mut model,
        ogeom::geom::SurfaceGeometry::Plane(ogeom::geom::PlaneSurface::new(
            ogeom::math::Plane::new(Frame::WORLD),
        )),
        // Lower arc first, then the upper against its sense: counter-
        // clockwise about +z, so the prism winds outward.
        &[vec![e_dn, e_up.reversed()]],
        T,
    )
    .unwrap()
    .shape;
    let solid = ogeom::algo::make_prism(
        &mut model,
        &face,
        ogeom::math::Vector::new(0.0, 0.0, 12.0),
        T,
    )
    .unwrap()
    .shape;
    let (original, recovered, import) = round_trip(model, solid);
    assert!(
        (recovered - original).abs() < original * 1e-6,
        "{recovered} against {original}"
    );
    assert_eq!(
        surface_kinds(import.document.model(), &import.solids[0])
            .iter()
            .filter(|k| **k == SurfaceKind::Extrusion)
            .count(),
        2,
        "both walls come back as extrusions"
    );
}

/// The surface kinds of a shape's faces.
fn surface_kinds(model: &Model, shape: &Shape) -> Vec<SurfaceKind> {
    ogeom::topo::explore_unique(model, shape, ogeom::topo::ShapeType::Face)
        .unwrap()
        .iter()
        .map(|f| {
            let Some(ogeom::topo::NodeData::Face(d)) = model.node(f).map(|n| n.data()) else {
                panic!("a face holds no face data");
            };
            ogeom::geom::Surface::kind(model.geometry().surface(d.surface).unwrap())
        })
        .collect()
}

/// A profile in the xz plane with a B-spline side, turned three quarters
/// of the way round the z axis: its curved wall a
/// surface of revolution of the spline, which goes out as one (120) over
/// the spline (126) with its start and terminate angles, and comes back
/// as one with the solid's volume.
#[test]
fn a_revolved_spline_profile_round_trips_through_its_surface_of_revolution() {
    let mut model = Model::new();
    let pts: Vec<Point> = (0..=8)
        .map(|i| {
            let t = f64::from(i) / 8.0;
            Point::new(4.0 + 1.5 * (core::f64::consts::PI * t).sin(), 0.0, 6.0 * t)
        })
        .collect();
    let side = ogeom::geom::Curve::from(
        ogeom::geom::fit::fit_points(&pts, 3, 1e-9, T)
            .unwrap()
            .curve,
    );
    let corners = [
        Point::new(4.0, 0.0, 0.0),
        Point::new(4.0, 0.0, 6.0),
        Point::new(1.0, 0.0, 6.0),
        Point::new(1.0, 0.0, 0.0),
    ];
    let v: Vec<Shape> = corners
        .iter()
        .map(|p| ogeom::algo::make_vertex(&mut model, *p).shape)
        .collect();
    let (s0, s1) = ogeom::geom::Curve3d::domain(&side);
    let wall = ogeom::algo::make_edge_between(&mut model, side, (s0, s1), &v[0], &v[1], T)
        .unwrap()
        .shape;
    let segment = |model: &mut Model, a: usize, b: usize| {
        let line = ogeom::geom::LineCurve::segment(corners[a], corners[b], T).unwrap();
        let range = ogeom::geom::Curve3d::domain(&line);
        ogeom::algo::make_edge_between(model, line.into(), range, &v[a], &v[b], T)
            .unwrap()
            .shape
    };
    let top = segment(&mut model, 1, 2);
    let inner = segment(&mut model, 2, 3);
    let bottom = segment(&mut model, 3, 0);
    let plane =
        ogeom::math::Plane::new(Frame::new(Point::ORIGIN, -Direction::Y, Direction::X, T).unwrap());
    let face = ogeom::algo::make_face_with_pcurves(
        &mut model,
        ogeom::geom::SurfaceGeometry::Plane(ogeom::geom::PlaneSurface::new(plane)),
        &[vec![wall, top, inner, bottom]],
        T,
    )
    .unwrap()
    .shape;
    let solid = ogeom::algo::make_revolution(
        &mut model,
        &face,
        ogeom::math::Axis {
            location: Point::ORIGIN,
            direction: Direction::Z,
        },
        1.5 * core::f64::consts::PI,
        T,
    )
    .unwrap()
    .shape;
    let before = surface_kinds(&model, &solid);
    assert!(before.contains(&SurfaceKind::Revolution), "{before:?}");
    let (original, recovered, import) = round_trip(model, solid);
    assert!(
        (recovered - original).abs() < original * 1e-6,
        "{recovered} against {original}"
    );
    let after = surface_kinds(import.document.model(), &import.solids[0]);
    assert_eq!(
        after
            .iter()
            .filter(|k| **k == SurfaceKind::Revolution)
            .count(),
        before
            .iter()
            .filter(|k| **k == SurfaceKind::Revolution)
            .count(),
        "{before:?} against {after:?}"
    );
}

#[test]
fn a_sphere_round_trips_through_its_seam_only_boundary() {
    // A full sphere's only written boundary is its seam, listed twice. The
    // poles have no curve to write. The reader recognises a seam-only loop
    // as the whole chart and rebuilds the natural face, degenerate boundary
    // and all. The two tessellations sample the chart differently, so the
    // comparison is against the closed form, at a chord fine enough that
    // the inscribed deficit stays inside the stated bound.
    let mut model = Model::new();
    let solid = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 7.0, T)
        .unwrap()
        .shape;
    let mut document = ogeom::doc::Document::over(std::mem::take(&mut model));
    document.add_part("ball", solid);
    let text = ogeom::io::write_iges(&document, T).unwrap();
    let import = ogeom::io::read_iges(&text, T).unwrap();
    let fine = Deflection::with_chord(0.01).unwrap();
    let recovered =
        ogeom::algo::volume_properties(import.document.model(), &import.solids[0], fine, T)
            .unwrap()
            .mass;
    let exact = 4.0 / 3.0 * core::f64::consts::PI * 7.0_f64.powi(3);
    // An inscribed tessellation at chord ε on radius r reads low by about
    // 3ε/r (every triangle sits up to ε inside), which at 0.01 on 7 is
    // 0.43%. The bound is that, with a third again for the sampling's
    // unevenness, and it is a statement about tessellation rather than
    // about the exchange: the same solid measured natively reads the same.
    let budget = 3.0 * 0.01 / 7.0 * (4.0 / 3.0);
    assert!(
        (recovered - exact).abs() < exact * budget,
        "{recovered} against {exact}, budget {budget}"
    );
}

#[test]
fn an_inch_file_scales_into_millimetres() {
    let mut model = Model::new();
    let solid = ogeom::algo::make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
        .unwrap()
        .shape;
    let mut document = ogeom::doc::Document::over(std::mem::take(&mut model));
    document.add_part("cube", solid);
    let text = ogeom::io::write_iges(&document, T).unwrap();
    // The same deck, restated in inches: the global units flag and name are
    // the only difference, and every coordinate must arrive 25.4 times as
    // large. The substitution is same-length so the fixed columns survive.
    let inch = text.replace(",2,2HMM,", ",1,2HIN,");
    assert_ne!(inch, text, "the substitution must strike");
    let import = ogeom::io::read_iges(&inch, T).unwrap();
    assert!((import.report.scale_mm - 25.4).abs() < 1e-12);
    let recovered = volume(import.document.model(), &import.solids[0]);
    let expected = 25.4_f64.powi(3);
    assert!(
        (recovered - expected).abs() < expected * 1e-9,
        "{recovered} against {expected}"
    );
}

#[test]
fn an_empty_deck_is_refused_by_name() {
    let err = ogeom::io::read_iges("not iges at all\n", T).unwrap_err();
    assert!(
        err.to_string().contains("section letter"),
        "unexpected message: {err}"
    );
}

/// A closed spline wall through exchange, both formats. The skinned
/// loft's wall is a clamped B-spline that closes on itself without being
/// periodic. Projections near the joining column land in either copy, and
/// only continuity chooses. The shared fitted-pcurve machinery unwraps by
/// *closure*, so what one reader learns the other has too, and the pin is
/// that both formats return the same solid, measured.
#[test]
fn f5_a_closed_spline_wall_survives_both_formats() {
    let mut model = Model::new();
    let profile = |model: &mut Model, z: f64, half: f64| -> Shape {
        let corners = [
            Point::new(-half, -half, z),
            Point::new(half, -half, z),
            Point::new(half, half, z),
            Point::new(-half, half, z),
        ];
        ogeom::algo::make_polygon(model, &corners, true, T)
            .unwrap()
            .shape
    };
    let a = profile(&mut model, 0.0, 8.0);
    let b = profile(&mut model, 5.0, 7.0);
    let c = profile(&mut model, 10.0, 5.0);
    let solid = ogeom::offset::make_loft_skinned(&mut model, &[a, b, c], 0.5, T)
        .unwrap()
        .shape;
    let original = volume(&model, &solid);
    let mut document = ogeom::doc::Document::over(std::mem::take(&mut model));
    document.add_part("loft", solid);

    let step = ogeom::io::read_step(&ogeom::io::write_step(&document, T).unwrap(), T).unwrap();
    let iges = ogeom::io::read_iges(&ogeom::io::write_iges(&document, T).unwrap(), T).unwrap();
    let via_step = volume(step.document.model(), &step.solids[0]);
    let via_iges = volume(iges.document.model(), &iges.solids[0]);

    // The fitted trims re-tessellate a hair differently from the native
    // face. The bound is that hair, and the two formats must agree with
    // each other far tighter than either must agree with the original.
    assert!(
        (via_step - original).abs() < original * 1e-3,
        "{via_step} against {original}"
    );
    assert!(
        (via_iges - original).abs() < original * 1e-3,
        "{via_iges} against {original}"
    );
    assert!(
        (via_step - via_iges).abs() < original * 1e-9,
        "one fix, two readers: {via_step} vs {via_iges}"
    );
}

/// Every wire of a recovered face closes in its own chart, both formats.
///
/// A periodic chart offers a branch per turn, all describing the same
/// points, and a reader that derives each edge's image on its own has no
/// reason to pick one over another. A drilled block whose bore wall's two
/// rims come back whole turns apart (the right circles, neither meeting
/// the seam the wire closes on) is caught by nothing downstream: the
/// tessellator folds into the chart and copes, and the error shows up only
/// as a tenth of a percent of volume.
///
/// Both readers slide each image by whole periods until its start meets
/// where the last one ended, and this walks the result to say so.
#[test]
fn a_recovered_wire_closes_in_its_own_chart() {
    use ogeom::geom::Curve2d as _;
    use ogeom::topo::{EdgeRepr, Orientation, ShapeType, explore_unique};

    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (30.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(15.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, frame, 4.0, 12.0, T)
        .unwrap()
        .shape;
    let solid = ogeom::boolean::cut(&mut model, &block, &drill, T)
        .unwrap()
        .shape;
    let mut document = ogeom::doc::Document::over(std::mem::take(&mut model));
    document.add_part("part", solid);

    let step = ogeom::io::read_step(&ogeom::io::write_step(&document, T).unwrap(), T).unwrap();
    let iges = ogeom::io::read_iges(&ogeom::io::write_iges(&document, T).unwrap(), T).unwrap();
    for (format, model, solid) in [
        ("step", step.document.model(), &step.solids[0]),
        ("iges", iges.document.model(), &iges.solids[0]),
    ] {
        let mut walls = 0;
        for face in explore_unique(model, solid, ShapeType::Face).unwrap() {
            let surface = model
                .node(&face)
                .and_then(|n| n.data().as_face())
                .map(|d| d.surface)
                .unwrap();
            if matches!(
                model.geometry().surface(surface),
                Some(ogeom::geom::SurfaceGeometry::Cylinder(_))
            ) {
                walls += 1;
            }
            for wire in model.ordered_children_of(&face).unwrap() {
                let mut first = None;
                let mut previous: Option<ogeom::math::Point2> = None;
                // Begun at an edge that is not the seam, so the first step
                // has somewhere to continue from: which column a seam's
                // occurrence takes is decided by the ring, as the
                // tessellator decides it, and not by any flag.
                let mut ring = model.ordered_children_of(&wire).unwrap();
                if let Some(at) = ring.iter().position(|e| {
                    let data = model.node(e).unwrap().data().as_edge().unwrap();
                    !matches!(
                        data.pcurve_for(surface, e.location()),
                        Some(EdgeRepr::Seam { .. })
                    )
                }) {
                    ring.rotate_left(at);
                }
                for edge in ring {
                    let back = edge.orientation() == Orientation::Reversed;
                    let data = model.node(&edge).unwrap().data().as_edge().unwrap();
                    let (image, range) = match data.pcurve_for(surface, edge.location()) {
                        Some(EdgeRepr::PCurve { curve, range, .. }) => {
                            (model.geometry().pcurve(*curve), *range)
                        }
                        // The column this occurrence takes: the one whose
                        // oriented start continues where the walk has got to.
                        Some(EdgeRepr::Seam {
                            forward,
                            reversed,
                            range,
                            ..
                        }) => {
                            let start = if back { range.1 } else { range.0 };
                            let nearest = |id| {
                                model
                                    .geometry()
                                    .pcurve(id)
                                    .and_then(|pc| pc.point_at(start, T).ok())
                                    .zip(previous)
                                    .map_or(f64::INFINITY, |(at, previous)| previous.distance(at))
                            };
                            let id = if nearest(*forward) <= nearest(*reversed) {
                                *forward
                            } else {
                                *reversed
                            };
                            (model.geometry().pcurve(id), *range)
                        }
                        _ => (None, (0.0, 0.0)),
                    };
                    let Some(image) = image else { continue };
                    let (start, end) = if back {
                        (range.1, range.0)
                    } else {
                        (range.0, range.1)
                    };
                    let at = image.point_at(start, T).unwrap();
                    let leaves = image.point_at(end, T).unwrap();
                    if let Some(previous) = previous {
                        assert!(
                            previous.distance(at) < 1e-6,
                            "{format}: a wire jumps {} across its chart",
                            previous.distance(at)
                        );
                    }
                    first.get_or_insert(at);
                    previous = Some(leaves);
                }
                if let (Some(first), Some(previous)) = (first, previous) {
                    assert!(
                        previous.distance(first) < 1e-6,
                        "{format}: a wire does not close in its chart, by {}",
                        previous.distance(first)
                    );
                }
            }
        }
        assert_eq!(walls, 1, "{format}: the bore's wall came back");
    }
}

/// A sheet on an offset of a bicubic bump, which has no analytic form, as
/// it stands and turned over: each goes out as an offset surface (140)
/// over the patch (128) and comes back as an offset surface with the
/// same area and facing the same way.
#[test]
fn an_offset_sheet_round_trips_as_an_offset_surface_either_way_up() {
    use ogeom::geom::{BSplineSurface, OffsetSurface, Surface as _, SurfaceGeometry};
    use ogeom::math::{ControlGrid, KnotVector};
    let mut points = Vec::new();
    for i in 0..4 {
        for j in 0..4 {
            let bump = if (1..3).contains(&i) && (1..3).contains(&j) {
                1.5
            } else {
                0.0
            };
            points.push(Point::new(f64::from(i) * 3.0, f64::from(j) * 3.0, bump));
        }
    }
    let patch = BSplineSurface::new(
        KnotVector::clamped_uniform(3, 4).unwrap(),
        KnotVector::clamped_uniform(3, 4).unwrap(),
        &ControlGrid::new(points, 4, 4).unwrap(),
        T,
    )
    .unwrap();
    let offset = OffsetSurface::new(SurfaceGeometry::BSpline(patch), 0.5).unwrap();
    assert!(offset.analytic(T).unwrap().is_none(), "the premise");
    let fine = Deflection {
        chord: 1e-4,
        ..Deflection::default()
    };
    for turned in [false, true] {
        let mut model = Model::new();
        let face = ogeom::algo::make_natural_face(
            &mut model,
            SurfaceGeometry::Offset(Box::new(offset.clone())),
        )
        .unwrap()
        .shape;
        let face = if turned { face.reversed() } else { face };
        let area = ogeom::algo::surface_properties(&model, &face, fine, T)
            .unwrap()
            .mass;
        let (_, facing) = ogeom::algo::face_normal(&model, &face, T).unwrap();
        let mut document = ogeom::doc::Document::over(model);
        document.add_part("sheet", face);
        let text = ogeom::io::write_iges(&document, T).unwrap();
        assert!(
            text.lines().any(|l| l.starts_with("     140")),
            "turned {turned}"
        );
        let import = ogeom::io::read_iges(&text, T).unwrap();
        let back = import.document.model();
        assert_eq!(import.sheets.len(), 1, "turned {turned}");
        let faces =
            ogeom::topo::explore_unique(back, &import.sheets[0], ogeom::topo::ShapeType::Face)
                .unwrap();
        assert_eq!(faces.len(), 1);
        let Some(ogeom::topo::NodeData::Face(d)) = back.node(&faces[0]).map(|n| n.data()) else {
            panic!("a face holds no face data");
        };
        assert_eq!(
            back.geometry().surface(d.surface).unwrap().kind(),
            SurfaceKind::Offset,
            "turned {turned}"
        );
        let again = ogeom::algo::surface_properties(back, &faces[0], fine, T)
            .unwrap()
            .mass;
        assert!(
            (again - area).abs() < area * 1e-9,
            "turned {turned}: {again} against {area}"
        );
        let (_, facing_back) = ogeom::algo::face_normal(back, &faces[0], T).unwrap();
        assert!(
            facing.z * facing_back.z > 0.0,
            "turned {turned}: {facing:?} against {facing_back:?}"
        );
    }
}

/// A part named past the eight columns of an entity label goes out with
/// its whole name as the solid's Name property (406 form 15) and comes
/// back under it; a name a string parameter cannot carry whole keeps to
/// the label.
#[test]
fn a_long_part_name_round_trips_through_its_name_property() {
    for (name, back) in [
        ("bracket_mount_left", "bracket_mount_left"),
        ("left, upper", "left, up"),
    ] {
        let mut model = Model::new();
        let solid = ogeom::algo::make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T)
            .unwrap()
            .shape;
        let mut document = ogeom::doc::Document::over(model);
        document.add_part(name, solid);
        let text = ogeom::io::write_iges(&document, T).unwrap();
        let import = ogeom::io::read_iges(&text, T).unwrap();
        let names: Vec<&str> = import
            .document
            .products()
            .map(|(_, p)| p.name.as_str())
            .collect();
        assert_eq!(names, vec![back], "{name}");
        assert!(
            import.report.skipped.is_empty(),
            "{name}: {:?}",
            import.report.skipped
        );
    }
}
