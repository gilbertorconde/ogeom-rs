//! §H of `docs/PLAN.md`: removing a set of faces, the wound closed from the
//! neighbours' own geometry, measured against the solids the features were
//! cut from, not against plausibility.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::SurfaceGeometry;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// Faces of `solid` whose surface satisfies the predicate.
fn faces_where(
    model: &Model,
    solid: &Shape,
    pred: impl Fn(&SurfaceGeometry) -> bool,
) -> Vec<Shape> {
    explore_unique(model, solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            model
                .node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(&pred)
        })
        .collect()
}

/// A through bore's wall removed: the rims are inner loops of the lid and
/// base, and the block comes back to the last bit of its exact volume:
/// no overshoot, because nothing is filled; the boundary is resewn.
#[test]
fn removing_a_bores_wall_makes_the_block_whole() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (30.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(15.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, frame, 4.0, 12.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &block, &drill, T)
        .unwrap()
        .shape;

    let walls = faces_where(&model, &drilled, |s| {
        matches!(s, SurfaceGeometry::Cylinder(_))
    });
    assert_eq!(walls.len(), 1, "one bore wall");
    let built = ogeom::boolean::remove_faces(&mut model, &drilled, &walls, T).unwrap();

    let healed = volume(&model, &built.shape);
    // Within float accumulation, not within a modelling tolerance: the lid
    // and base come from the boolean's rebuild and carry its arithmetic. A
    // wrong closure would miss by the bore's volume, eight hundred times the
    // bound.
    assert!((healed - 6000.0).abs() < 1e-6, "the block, whole: {healed}");
    assert!(built.history.is_deleted(&walls[0]));
}

/// A chamfer band removed: the band interrupts the top and side faces, whose
/// planes re-intersect in the edge the chamfer replaced, and the end caps'
/// edges extend to the recovered corners.
#[test]
fn removing_a_chamfer_restores_the_sharp_box() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 15.0, 10.0), T)
        .unwrap()
        .shape;
    let edge = edge_near(&model, &block, Point::new(10.0, 0.0, 10.0));
    let chamfered = ogeom::fillet::chamfer_edge(&mut model, &block, &edge, 3.0, T)
        .unwrap()
        .shape;
    let before = volume(&model, &chamfered);
    assert!(before < 3000.0, "the chamfer took material: {before}");

    // The chamfer band is the one plane that is neither axis-aligned wall
    // nor lid: its normal has both y and z.
    let bands = faces_where(&model, &chamfered, |s| {
        if let SurfaceGeometry::Plane(p) = s {
            let n = p.plane().frame().z().vector();
            n.y.abs() > 0.1 && n.z.abs() > 0.1
        } else {
            false
        }
    });
    assert_eq!(bands.len(), 1, "one chamfer band");
    let built = ogeom::boolean::remove_faces(&mut model, &chamfered, &bands, T).unwrap();

    let healed = volume(&model, &built.shape);
    assert!(
        (healed - 3000.0).abs() < 1e-9,
        "the sharp box, exactly: {healed}"
    );
    // And it is a box again: six faces, twelve edges, eight vertices.
    let counts = |kind| explore_unique(&model, &built.shape, kind).unwrap().len();
    assert_eq!(counts(ShapeType::Face), 6);
    assert_eq!(counts(ShapeType::Edge), 12);
    assert_eq!(counts(ShapeType::Vertex), 8);
}

/// A fillet band removed: the same wound with a cylindrical band, and the
/// caps' arcs replaced by their own straight edges extended to the corner.
#[test]
fn removing_a_fillet_restores_the_sharp_box() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 15.0, 10.0), T)
        .unwrap()
        .shape;
    let edge = edge_near(&model, &block, Point::new(10.0, 0.0, 10.0));
    let filleted = ogeom::fillet::fillet_edge(&mut model, &block, &edge, 3.0, T)
        .unwrap()
        .shape;

    let bands = faces_where(&model, &filleted, |s| {
        matches!(s, SurfaceGeometry::Cylinder(_))
    });
    assert_eq!(bands.len(), 1, "one fillet band");
    let built = ogeom::boolean::remove_faces(&mut model, &filleted, &bands, T).unwrap();

    let healed = volume(&model, &built.shape);
    assert!(
        (healed - 3000.0).abs() < 1e-9,
        "the sharp box, exactly: {healed}"
    );
}

/// The refusals, by name.
#[test]
fn impossible_removals_are_refused_by_name() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let faces = explore_unique(&model, &block, ShapeType::Face).unwrap();

    let err = ogeom::boolean::remove_faces(&mut model, &block, &[], T).unwrap_err();
    assert!(err.to_string().contains("nothing to remove"), "{err}");

    let err = ogeom::boolean::remove_faces(&mut model, &block, &faces, T).unwrap_err();
    assert!(err.to_string().contains("nothing remains"), "{err}");

    // Removing one wall of a plain box: the wound has no feature geometry to
    // close it: the neighbours meet at right angles already, and nothing
    // stands in for the missing face.
    let err = ogeom::boolean::remove_faces(&mut model, &block, &faces[..1], T).unwrap_err();
    assert!(
        !err.to_string().is_empty(),
        "a wall removal fails with a reason, not silently"
    );
}

/// The box edge whose midpoint is nearest `at`.
fn edge_near(model: &Model, solid: &Shape, at: Point) -> Shape {
    use ogeom::geom::Curve3d as _;
    let mut best: Option<(f64, Shape)> = None;
    for edge in explore_unique(model, solid, ShapeType::Edge).unwrap() {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            continue;
        };
        let p = geometry
            .point_at(f64::midpoint(range.0, range.1), T)
            .unwrap();
        let d = p.distance(at);
        if best.as_ref().is_none_or(|(held, _)| d < *held) {
            best = Some((d, edge));
        }
    }
    best.unwrap().1
}

/// Two separate features named in one call remove as two wounds.
///
/// Two fillets on opposite edges share nothing; classifying their ring
/// edges together used to elect "the sides" across both features and die
/// recovering a nonsense edge. Grouped by shared edges, each feature runs
/// the whole machinery on the previous result, and the box comes back
/// sharp: exactly, both wounds.
#[test]
fn two_separate_fillets_remove_in_one_call() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let mut solid = block;
    for target in [Point::new(10.0, 5.0, 10.0), Point::new(0.0, 5.0, 0.0)] {
        let edge = edge_near(&model, &solid, target);
        solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, 2.0, T)
            .unwrap()
            .shape;
    }
    // The two blend faces: the cylindrical ones.
    let mut blends = Vec::new();
    for face in explore_unique(&model, &solid, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        let surface = model.geometry().surface(data.surface).unwrap();
        if matches!(surface, ogeom::geom::SurfaceGeometry::Cylinder(_)) {
            blends.push(face);
        }
    }
    assert_eq!(blends.len(), 2, "two fillets leave two cylindrical faces");
    let restored = ogeom::boolean::remove_faces(&mut model, &solid, &blends, T)
        .unwrap()
        .shape;
    assert!(
        ogeom::algo::check(&model, &restored, T).unwrap().is_valid(),
        "the restored box is a valid solid"
    );
    let volume =
        ogeom::algo::volume_properties(&model, &restored, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    assert!(
        (volume - 1000.0).abs() < 1e-6,
        "both wounds close exactly: {volume} against 1000"
    );
}

/// Two blends meeting at a box corner (`fillet_edges` on two top edges,
/// the later band trimmed against the earlier) removed in one call: each
/// band recovers its own crease from its side planes, the two creases meet
/// where one pierces the other's side, that corner is one vertex for both,
/// and the box comes back to the last bit.
#[test]
fn two_meeting_fillets_remove_in_one_call() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let a = edge_near(&model, &block, Point::new(10.0, 20.0, 10.0));
    let b = edge_near(&model, &block, Point::new(20.0, 10.0, 10.0));
    let met = ogeom::fillet::fillet_edges(&mut model, &block, &[a, b], 2.0, T)
        .unwrap()
        .shape;
    let bands: Vec<Shape> = explore_unique(&model, &met, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::Cylinder(_))
            )
        })
        .collect();
    assert_eq!(bands.len(), 2);
    let restored = ogeom::boolean::remove_faces(&mut model, &met, &bands, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &restored, T).unwrap().is_valid());
    assert_eq!(
        explore_unique(&model, &restored, ShapeType::Face)
            .unwrap()
            .len(),
        6
    );
    assert!((volume(&model, &restored) - 4000.0).abs() < 1e-6);
}

/// The same corner blended one edge at a time: the second wedge's cap
/// stands flush against the first band, and the feature is the two bands
/// *and* that cap. The cap borders one wall and two removed faces; it
/// joins its band's crease, whose end is then where the crease pierces the
/// other band's side wall: the shared corner.
#[test]
fn two_flush_fillets_and_their_cap_remove_in_one_call() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let a = edge_near(&model, &block, Point::new(10.0, 20.0, 10.0));
    let first = ogeom::fillet::fillet_edge(&mut model, &block, &a, 2.0, T)
        .unwrap()
        .shape;
    let b = edge_near(&model, &first, Point::new(20.0, 10.0, 10.0));
    let flush = ogeom::fillet::fillet_edge(&mut model, &first, &b, 2.0, T)
        .unwrap()
        .shape;
    let faces = explore_unique(&model, &flush, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 9, "six walls, two bands, one cap");
    let feature: Vec<Shape> = faces
        .into_iter()
        .filter(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::Cylinder(_))
            ) || explore_unique(&model, f, ShapeType::Edge).unwrap().len() == 3
        })
        .collect();
    assert_eq!(feature.len(), 3);
    let restored = ogeom::boolean::remove_faces(&mut model, &flush, &feature, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &restored, T).unwrap().is_valid());
    assert_eq!(
        explore_unique(&model, &restored, ShapeType::Face)
            .unwrap()
            .len(),
        6
    );
    assert!((volume(&model, &restored) - 4000.0).abs() < 1e-6);
}

/// The issue's own acceptance: two chamfers meeting at a box corner,
/// removed in one call. One at a time they stand flush, the second's
/// triangular cap against the first's plane; the feature is both chamfer
/// planes and the cap.
#[test]
fn two_chamfers_meeting_at_a_corner_remove_in_one_call() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let a = edge_near(&model, &block, Point::new(10.0, 20.0, 10.0));
    let first = ogeom::fillet::chamfer_edge(&mut model, &block, &a, 2.0, T)
        .unwrap()
        .shape;
    let b = edge_near(&model, &first, Point::new(20.0, 10.0, 10.0));
    let both = ogeom::fillet::chamfer_edge(&mut model, &first, &b, 2.0, T)
        .unwrap()
        .shape;
    let faces = explore_unique(&model, &both, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 9, "six walls, two chamfer planes, one cap");
    let feature: Vec<Shape> = faces
        .into_iter()
        .filter(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            let Some(SurfaceGeometry::Plane(p)) = model.geometry().surface(data.surface) else {
                return true;
            };
            let n = p.plane().normal().vector();
            let off_axis = !(n.x.abs() > 0.999 || n.y.abs() > 0.999 || n.z.abs() > 0.999);
            off_axis || explore_unique(&model, f, ShapeType::Edge).unwrap().len() == 3
        })
        .collect();
    assert_eq!(feature.len(), 3);
    let restored = ogeom::boolean::remove_faces(&mut model, &both, &feature, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &restored, T).unwrap().is_valid());
    assert_eq!(
        explore_unique(&model, &restored, ShapeType::Face)
            .unwrap()
            .len(),
        6
    );
    assert!((volume(&model, &restored) - 4000.0).abs() < 1e-6);
}

/// A rim blend's wound is the whole ring it took out of each neighbour,
/// and it closes on the circle those neighbours meet along.
///
/// Three rims, each blended and then taken away again: a drum's top, where
/// the cap's *whole outer boundary* was the wound and the cap must grow
/// back; a bore's mouth in a plate, where the wound is an inner ring of
/// the top and the whole bottom of the bore's wall; and a boss's seat,
/// where the blend is additive and the wound sits between the plate's top
/// and the post. Each comes back to the solid it was cut from, face for
/// face and to the volume the mesh can measure.
///
/// Each with a fillet and with a chamfer, since the wound is the same
/// whether the band is a torus or a cone.
///
/// This is the wound an earlier note called "a neighbour meeting itself".
/// A whole ring taken out of a neighbour is closed one of two ways, and
/// only the neighbours' surfaces say which: a bore's two mouths sit in
/// faces that never meet, so the rings are dropped and the faces grow
/// over them, while a rim blend's neighbours meet along the very circle
/// it replaced. The wall's own seam then reaches that circle, which is
/// where the circle is cut, and the seam extends to meet it.
#[test]
fn a_rim_blend_removes_and_its_rim_comes_back() {
    for (case, bevel) in [
        ("drum", false),
        ("drum", true),
        ("mouth", false),
        ("mouth", true),
        ("boss", false),
        ("boss", true),
    ] {
        let mut model = Model::new();
        let r = 1.0;
        let (sharp, rim_at) = match case {
            "drum" => {
                let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
                    .unwrap()
                    .shape;
                (drum, Point::new(0.0, 5.0, 10.0))
            }
            "mouth" => {
                let block = ogeom::algo::make_box(
                    &mut model,
                    Frame::new(Point::new(-10.0, -10.0, 0.0), Direction::Z, Direction::X, T)
                        .unwrap(),
                    (20.0, 20.0, 10.0),
                    T,
                )
                .unwrap()
                .shape;
                let seat =
                    Frame::new(Point::new(0.0, 0.0, -1.0), Direction::Z, Direction::X, T).unwrap();
                let drill = ogeom::algo::make_cylinder(&mut model, seat, 4.0, 12.0, T)
                    .unwrap()
                    .shape;
                let bored = ogeom::boolean::cut(&mut model, &block, &drill, T)
                    .unwrap()
                    .shape;
                (bored, Point::new(0.0, 4.0, 10.0))
            }
            _ => {
                let plate = ogeom::algo::make_box(
                    &mut model,
                    Frame::new(Point::new(-10.0, -10.0, 0.0), Direction::Z, Direction::X, T)
                        .unwrap(),
                    (20.0, 20.0, 5.0),
                    T,
                )
                .unwrap()
                .shape;
                let seat =
                    Frame::new(Point::new(0.0, 0.0, 5.0), Direction::Z, Direction::X, T).unwrap();
                let post = ogeom::algo::make_cylinder(&mut model, seat, 4.0, 6.0, T)
                    .unwrap()
                    .shape;
                let boss = ogeom::boolean::fuse(&mut model, &plate, &post, T)
                    .unwrap()
                    .shape;
                (boss, Point::new(0.0, 4.0, 5.0))
            }
        };
        // Measured off the mesh at a fine chord, because the blend's own
        // boolean leaves the faces it did not touch split at their seams,
        // and a solid whose discs are two arcs apiece no longer takes the
        // exact integrator's path. The reach quoted below is the mesh's,
        // and it is two orders finer than the feature being removed.
        let fine = Deflection::with_chord(1e-4).unwrap();
        let measure = |model: &Model, shape: &Shape| {
            ogeom::algo::volume_properties(model, shape, fine, T)
                .unwrap()
                .mass
        };
        let was = measure(&model, &sharp);
        let faces_before = explore_unique(&model, &sharp, ShapeType::Face)
            .unwrap()
            .len();

        let rim = edge_near(&model, &sharp, rim_at);
        let made = if bevel {
            ogeom::fillet::chamfer_edge(&mut model, &sharp, &rim, r, T)
        } else {
            ogeom::fillet::fillet_edge(&mut model, &sharp, &rim, r, T)
        };
        let blended = made
            .unwrap_or_else(|e| panic!("{case} bevel {bevel}: {e}"))
            .shape;
        let band = faces_where(&model, &blended, |s| {
            if bevel {
                matches!(s, SurfaceGeometry::Cone(_))
            } else {
                matches!(s, SurfaceGeometry::Torus(_))
            }
        });
        assert_eq!(band.len(), 1, "{case} bevel {bevel}: one band");

        let back = ogeom::boolean::remove_faces(&mut model, &blended, &band, T)
            .unwrap_or_else(|e| panic!("{case} bevel {bevel}: {e}"))
            .shape;
        let diagnosis = ogeom::algo::check(&model, &back, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "{case} bevel {bevel}: {:?}",
            diagnosis.problems
        );
        assert_eq!(
            explore_unique(&model, &back, ShapeType::Face)
                .unwrap()
                .len(),
            faces_before,
            "{case} bevel {bevel}: the faces the feature interrupted are whole again"
        );
        let now = measure(&model, &back);
        assert!(
            (now - was).abs() < was * 1e-4,
            "{case} bevel {bevel}: {now} against the solid it was cut from, {was}"
        );
    }
}

/// A stadium prism: a rectangle with semicircular ends, extruded up.
fn stadium(model: &mut Model, length: f64, r: f64, height: f64) -> Shape {
    use ogeom::geom::{CircleCurve, Curve, Curve3d as _, LineCurve, PlaneSurface};
    use ogeom::math::{Circle, Plane, Vector};
    let tl = ogeom::algo::make_vertex(model, Point::new(0.0, r, 0.0)).shape;
    let tr = ogeom::algo::make_vertex(model, Point::new(length, r, 0.0)).shape;
    let br = ogeom::algo::make_vertex(model, Point::new(length, -r, 0.0)).shape;
    let bl = ogeom::algo::make_vertex(model, Point::new(0.0, -r, 0.0)).shape;
    let arc = |model: &mut Model, centre: Point, x: Direction, from: &Shape, to: &Shape| {
        let frame = Frame::new(centre, Direction::Z, x, T).unwrap();
        let curve = Curve::Circle(CircleCurve::new(Circle::new(frame, r, T).unwrap()));
        ogeom::algo::make_edge_between(model, curve, (0.0, core::f64::consts::PI), from, to, T)
            .unwrap()
            .shape
    };
    let seg = |model: &mut Model, from: (&Shape, Point), to: (&Shape, Point)| {
        let curve = Curve::Line(LineCurve::segment(from.1, to.1, T).unwrap());
        let domain = curve.domain();
        ogeom::algo::make_edge_between(model, curve, domain, from.0, to.0, T)
            .unwrap()
            .shape
    };
    let top = seg(
        model,
        (&tl, Point::new(0.0, r, 0.0)),
        (&tr, Point::new(length, r, 0.0)),
    );
    let right = arc(
        model,
        Point::new(length, 0.0, 0.0),
        Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
        &br,
        &tr,
    );
    let bottom = seg(
        model,
        (&br, Point::new(length, -r, 0.0)),
        (&bl, Point::new(0.0, -r, 0.0)),
    );
    let left = arc(model, Point::new(0.0, 0.0, 0.0), Direction::Y, &tl, &bl);
    let plane = PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom::algo::make_face_with_pcurves(
        model,
        plane.into(),
        &[vec![left, bottom.reversed(), right, top.reversed()]],
        T,
    )
    .unwrap()
    .shape;
    ogeom::algo::make_prism(model, &face, Vector::new(0.0, 0.0, height), T)
        .unwrap()
        .shape
}

/// The edges of `solid` whose every vertex stands at `height`.
fn edges_at_height(model: &Model, solid: &Shape, height: f64) -> Vec<Shape> {
    use ogeom::topo::{Filter, explore};
    explore_unique(model, solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            explore(model, e, Filter::OfType(ShapeType::Vertex))
                .unwrap()
                .iter()
                .all(|v| {
                    model
                        .node(v)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .zip(v.transform(model.datums()).ok())
                        .is_some_and(|(p, placed)| (placed.apply(p).z - height).abs() < 1e-9)
                })
        })
        .collect()
}

/// A tangent chain of blends (a stadium's whole top rim rounded in one
/// call, two straight bands and two semicircular ones meeting flush)
/// removed in one call. At a tangent junction neither band's crease
/// pierces the other's wall: the straight crease grazes the round wall,
/// and the round crease grazes the flat one. The corner is where the two
/// creases touch, and the cross-section edge the two bands share says
/// exactly where: its foot on either crease. The round creases are
/// circles, and a band standing across a circle's seam is read as one run
/// about its own centre rather than its complement, so each end wall grows
/// back its outer half and not its inner.
#[test]
fn a_tangent_chain_of_blends_removes_in_one_call() {
    let mut model = Model::new();
    let (length, r, height) = (10.0, 5.0, 4.0);
    let solid = stadium(&mut model, length, r, height);
    let before = volume(&model, &solid);
    let rim = edges_at_height(&model, &solid, height);
    assert_eq!(rim.len(), 4, "the stadium's top rim has four edges");
    let blended = ogeom::fillet::fillet_edges(&mut model, &solid, &rim, 1.0, T)
        .unwrap()
        .shape;
    let bands = faces_where(&model, &blended, |s| match s {
        SurfaceGeometry::Torus(_) => true,
        SurfaceGeometry::Cylinder(c) => (c.cylinder().radius() - 1.0).abs() < 1e-9,
        _ => false,
    });
    assert_eq!(bands.len(), 4, "two straight bands and two round ones");
    let restored = ogeom::boolean::remove_faces(&mut model, &blended, &bands, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &restored, T).unwrap().is_valid());
    assert_eq!(
        explore_unique(&model, &restored, ShapeType::Face)
            .unwrap()
            .len(),
        6,
        "the top, the bottom and four walls"
    );
    let shell = explore_unique(&model, &restored, ShapeType::Shell)
        .unwrap()
        .remove(0);
    assert!(ogeom::algo::is_shell_closed(&model, &shell).unwrap());
    assert!((volume(&model, &restored) - before).abs() < 1e-6);
}

/// The whole rim of a box top rounded in one call (four bands meeting at
/// four corners) and removed in one call: every corner is where one
/// crease pierces the wall of the next, and the box comes back sharp.
#[test]
fn a_loop_of_four_fillets_removes_in_one_call() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let rim = edges_at_height(&model, &block, 10.0);
    assert_eq!(rim.len(), 4);
    let blended = ogeom::fillet::fillet_edges(&mut model, &block, &rim, 2.0, T)
        .unwrap()
        .shape;
    let bands = faces_where(&model, &blended, |s| {
        matches!(s, SurfaceGeometry::Cylinder(_))
    });
    assert_eq!(bands.len(), 4);
    let restored = ogeom::boolean::remove_faces(&mut model, &blended, &bands, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &restored, T).unwrap().is_valid());
    assert_eq!(
        explore_unique(&model, &restored, ShapeType::Face)
            .unwrap()
            .len(),
        6
    );
    assert!((volume(&model, &restored) - 4000.0).abs() < 1e-6);
}
