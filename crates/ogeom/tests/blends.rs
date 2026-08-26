//! §10's tail: what a blend achieved, measured; blends between faces that
//! share no edge; edges whose envelope has no closed form; and the corner
//! where three of them meet.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The edge of `shape` whose midpoint is nearest `near`.
fn edge_near(model: &Model, shape: &Shape, near: Point) -> Shape {
    use ogeom::geom::Curve3d as _;
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| {
            let mid = |e: &Shape| {
                let data = model.node(e).unwrap().data().as_edge().unwrap();
                let ogeom::topo::EdgeRepr::Curve3d { curve, range, .. } = data.curve3d().unwrap()
                else {
                    unreachable!()
                };
                model
                    .geometry()
                    .curve(*curve)
                    .unwrap()
                    .point_at(f64::midpoint(range.0, range.1), T)
                    .unwrap()
                    .distance(near)
            };
            mid(a)
                .partial_cmp(&mid(b))
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .expect("some edge")
}

/// The planar face of `shape` whose plane passes through `on` and whose own
/// vertices bracket it.
fn planar_face_at(model: &Model, shape: &Shape, on: Point) -> Shape {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let ogeom::topo::NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            let Some(ogeom::geom::SurfaceGeometry::Plane(plane)) =
                model.geometry().surface(data.surface)
            else {
                return false;
            };
            if plane.plane().distance_to(on).abs() > 1e-9 {
                return false;
            }
            let mut bound = ogeom::math::Aabb::EMPTY;
            for v in explore_unique(model, f, ShapeType::Vertex).unwrap() {
                bound = bound.with_point(model.node(&v).unwrap().data().as_vertex().unwrap().point);
            }
            bound.expanded(1e-6).contains(on)
        })
        .expect("a planar face there")
}

#[test]
fn a_fillet_reports_its_own_tangency_instead_of_claiming_it() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 12.0), T)
        .unwrap()
        .shape;
    let edge = edge_near(&model, &block, Point::new(20.0, 0.0, 12.0));
    let blended = ogeom::fillet::fillet_edge(&mut model, &block, &edge, 2.0, T)
        .unwrap()
        .shape;

    // The blend is the one cylindrical face on the result.
    let blend = explore_unique(&model, &blended, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let ogeom::topo::NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            matches!(
                model.geometry().surface(data.surface),
                Some(ogeom::geom::SurfaceGeometry::Cylinder(_))
            )
        })
        .expect("the rolling ball left a cylinder");

    let contacts = ogeom::fillet::analyse_blend(&model, &blended, &blend, 9, T).unwrap();
    assert_eq!(contacts.len(), 4, "two tangency edges and two end caps");
    // The two long edges are the tangency lines: smooth to rounding. The
    // two ends are the cap arcs, where the blend meets a face it is *not*
    // tangent to — a right angle, and it should say so.
    let mut smooth = 0;
    let mut square = 0;
    for contact in &contacts {
        assert!(
            contact.gap < 1e-9,
            "the shared edge lies on both surfaces: {}",
            contact.gap
        );
        if contact.tangency_error < 1e-9 {
            smooth += 1;
        } else if (contact.tangency_error - core::f64::consts::FRAC_PI_2).abs() < 1e-9 {
            square += 1;
        }
    }
    assert_eq!(
        (smooth, square),
        (2, 2),
        "two tangent joins, two square ones: {contacts:?}"
    );
}

#[test]
fn a_blend_bridges_two_faces_that_share_no_edge() {
    // A step: a tall block and a low one side by side, their vertical wall
    // and horizontal lid meeting at no edge at all. The rolling ball still
    // has a seat — it touches both — and the blend is the fillet that seat
    // implies.
    let mut model = Model::new();
    let tall = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let low = ogeom::algo::make_box(
        &mut model,
        Frame::new(Point::new(10.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
        (20.0, 20.0, 10.0),
        T,
    )
    .unwrap()
    .shape;
    let step = ogeom::boolean::fuse(&mut model, &tall, &low, T)
        .unwrap()
        .shape;

    let wall = planar_face_at(&model, &step, Point::new(10.0, 10.0, 15.0));
    let lid = planar_face_at(&model, &step, Point::new(20.0, 10.0, 10.0));

    let blended = ogeom::fillet::blend_faces(&mut model, &step, &wall, &lid, 4.0, T)
        .unwrap()
        .shape;
    let volume =
        ogeom::algo::volume_properties(&model, &blended, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    // The step is 10*20*20 + 20*20*10 = 8000, and its inner corner is
    // concave: the ball rolls in the notch, so the blend *fills* it with
    // what a square corner would have held minus the quarter disc,
    // (r^2 - pi r^2 / 4), along the 20 of run.
    let r: f64 = 4.0;
    let filled = r.mul_add(r, -(core::f64::consts::PI * r * r / 4.0)) * 20.0;
    assert!(
        (volume - (8000.0 + filled)).abs() < 8000.0 * 2e-3,
        "the notch is filled, not cut: {volume} against {}",
        8000.0 + filled
    );
}

/// B2 — the corner where three blends meet. Three edges of a box are
/// filleted in sequence at one vertex, and the leftover spike is rounded by
/// the A5 tool: the corner block less the ball. The result is measured
/// against a closed form derived independently, by inclusion–exclusion over
/// the corner cube: within the cube every fillet prism's removal lies inside
/// the spike's, so the removed volume is three prism runs *outside* the cube
/// plus the spike itself —
///
///   V = 10³ − 3(1 − π/4) r² (10 − r) − r³ + πr³/6
///
/// which for r = 3 is 784 + 51.75π. The blend is tangent to everything it
/// rounds by construction — each contact a chart-degenerate curve or a
/// vertex of the tool's own patch — and this test is the corner family's
/// pin: it exercises A6, the tangential set-aside, the degeneracy splits,
/// and the tolerance-carrying welds at once.
#[test]
fn b2_three_fillets_and_the_corner_tool_round_the_vertex() {
    let mut model = Model::new();
    let r = 3.0;
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let mut solid = block;
    for target in [
        Point::new(10.0, 10.0, 5.0),
        Point::new(10.0, 5.0, 10.0),
        Point::new(5.0, 10.0, 10.0),
    ] {
        let edge = edge_near(&model, &solid, target);
        solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, r, T)
            .unwrap()
            .shape;
    }

    let at = |p: Point| Frame::new(p, Direction::Z, Direction::X, T).unwrap();
    let corner = Point::new(10.0 - r, 10.0 - r, 10.0 - r);
    let cblock = ogeom::algo::make_box(&mut model, at(corner), (r, r, r), T)
        .unwrap()
        .shape;
    let ball = ogeom::algo::make_sphere(&mut model, at(corner), r, T)
        .unwrap()
        .shape;
    let tool = ogeom::boolean::cut(&mut model, &cblock, &ball, T)
        .unwrap()
        .shape;
    let rounded = ogeom::boolean::cut(&mut model, &solid, &tool, T)
        .unwrap()
        .shape;

    assert!(
        ogeom::algo::check(&model, &rounded, T).unwrap().is_valid(),
        "the rounded corner is a valid solid"
    );
    let pi = core::f64::consts::PI;
    let want =
        1000.0 - 3.0 * (1.0 - pi / 4.0) * r * r * (10.0 - r) - r * r * r + pi * r * r * r / 6.0;
    let mut previous = f64::INFINITY;
    for chord in [1e-3, 1e-4] {
        let fine = ogeom::mesh::Deflection::with_chord(chord).unwrap();
        let measured = ogeom::algo::volume_properties(&model, &rounded, fine, T)
            .unwrap()
            .mass;
        let error = (measured - want).abs() / want;
        assert!(
            error < previous,
            "refining the mesh brings the measurement closer: {measured} vs {want}"
        );
        // The curved area is three band runs and the octant; the inscribed
        // deficit at chord δ runs to a few δ/r of the curved volume share.
        assert!(
            error < chord * 2.0,
            "the vertex blend against its closed form at chord {chord}: \
             {measured} vs {want}"
        );
        previous = error;
    }
}

/// The promoted corner tool: `round_vertex` reproduces the B2 closed form.
///
/// Same three fillets, same corner, same inclusion–exclusion reference —
/// but the ball-and-block construction now lives in the fillet crate with
/// its own refusals, instead of being spelled out per call site.
#[test]
fn round_vertex_reproduces_the_b2_closed_form() {
    let mut model = Model::new();
    let r = 3.0;
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    // The corner is captured while it still exists: the fillets consume the
    // tip, and the promoted tool reads the corner's planes from wherever the
    // vertex's *point* says they are.
    let vertex = vertex_near(&model, &block, Point::new(10.0, 10.0, 10.0));
    let mut solid = block;
    for target in [
        Point::new(10.0, 10.0, 5.0),
        Point::new(10.0, 5.0, 10.0),
        Point::new(5.0, 10.0, 10.0),
    ] {
        let edge = edge_near(&model, &solid, target);
        solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, r, T)
            .unwrap()
            .shape;
    }
    let rounded = ogeom::fillet::round_vertex(&mut model, &solid, &vertex, r, T)
        .unwrap()
        .shape;
    assert!(
        ogeom::algo::check(&model, &rounded, T).unwrap().is_valid(),
        "the rounded corner is a valid solid"
    );
    let expected = 784.0 + 51.75 * core::f64::consts::PI;
    for chord in [1e-3, 1e-4] {
        let fine = ogeom::mesh::Deflection::with_chord(chord).unwrap();
        let measured = ogeom::algo::volume_properties(&model, &rounded, fine, T)
            .unwrap()
            .mass;
        let error = (measured - expected).abs() / expected;
        assert!(
            error < chord * 2.0,
            "round_vertex against the closed form at chord {chord}: \
             {measured} vs {expected} ({error:.2e})"
        );
    }
}

/// The refusals name their families: a curved-edged corner and an oblique
/// one both belong to the setback construction, and say so.
#[test]
fn round_vertex_refuses_the_setback_family_by_name() {
    let mut model = Model::new();
    // A cylinder's rim vertex has a curved edge: refused as curved.
    let cyl = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
        .unwrap()
        .shape;
    let v = vertex_near(&model, &cyl, Point::new(5.0, 0.0, 10.0));
    let err = ogeom::fillet::round_vertex(&mut model, &cyl, &v, 1.0, T)
        .unwrap_err()
        .to_string();
    assert!(
        err.contains("still owed") || err.contains("exactly three"),
        "the curved corner names its family: {err}"
    );
}

/// The recon shape of issue #18 step 2: a box grooved by a tilted drum,
/// whose creases are ellipse arcs cut open by the box sides and split
/// again by the cylinder's own seam.
fn grooved_block(model: &mut Model) -> Shape {
    use ogeom::math::Vector;
    let block = ogeom::algo::make_box(model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let tilt = 0.35_f64;
    let axis = Direction::new(Vector::new(0.0, tilt.cos(), -tilt.sin()), T).unwrap();
    let frame = Frame::new(Point::new(10.0, -5.0, 11.5), axis, Direction::X, T).unwrap();
    let drum = ogeom::algo::make_cylinder(model, frame, 4.0, 30.0, T)
        .unwrap()
        .shape;
    ogeom::boolean::cut(model, &block, &drum, T).unwrap().shape
}

#[test]
fn an_open_seat_ends_in_run_out_caps() {
    // The bottom crease of the grooved block is an ellipse arc that meets
    // the box wall at both ends: an open seat. The marched band is trimmed
    // to the edge's own window and capped in each end section's own plane —
    // material comes off, the blend rides both hosts tangentially, and the
    // caps stand square to the band, which is what a run-out is.
    let mut model = Model::new();
    let grooved = grooved_block(&mut model);
    let before =
        ogeom::algo::volume_properties(&model, &grooved, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    let arc = edge_near(&model, &grooved, Point::new(10.0, 14.84, 0.0));
    let built = ogeom::fillet::fillet_edge(&mut model, &grooved, &arc, 1.0, T).unwrap();
    let after =
        ogeom::algo::volume_properties(&model, &built.shape, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    let removed = before - after;
    assert!(
        removed > 1.0 && removed < before * 0.05,
        "a run-out fillet removes a sliver, not a bite: {removed}"
    );

    // The blend face is the fitted band; its rails ride the hosts
    // tangentially and its end arcs stand square in the caps.
    use ogeom::topo::NodeData;
    let blend = explore_unique(&model, &built.shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let NodeData::Face(d) = model.node(f).unwrap().data() else {
                return false;
            };
            matches!(
                model.geometry().surface(d.surface),
                Some(ogeom::geom::SurfaceGeometry::BSpline(_))
            )
        })
        .expect("the fitted band is a face of the result");
    let contacts = ogeom::fillet::analyse_blend(&model, &built.shape, &blend, 15, T).unwrap();
    let mut smooth = 0;
    let mut square = 0;
    for c in &contacts {
        assert!(c.gap < 1e-3, "a contact stands off its edge: {}", c.gap);
        if c.tangency_error < 5e-3 {
            smooth += 1;
        } else if (c.tangency_error - core::f64::consts::FRAC_PI_2).abs() < 1e-3 {
            square += 1;
        }
    }
    assert_eq!(
        (smooth, square),
        (2, 2),
        "two tangent rails and two square caps: {contacts:?}"
    );
}

#[test]
fn a_seam_split_crease_arc_rounds_with_run_out_caps() {
    // The top crease is split by the cylinder's own seam into two arcs
    // sharing a mid-scoop vertex. The seat probe used to die on these —
    // the reconstructed loop's midpoint stands in cut-away territory —
    // before the march could speak. Probed and seated on the crease
    // itself, the arc marches its seat and lands as a capped blend.
    let mut model = Model::new();
    let grooved = grooved_block(&mut model);
    let arc = edge_near(&model, &grooved, Point::new(7.3, 7.7, 10.0));
    let built = ogeom::fillet::fillet_edge(&mut model, &grooved, &arc, 1.0, T).unwrap();
    // The result still meshes as one closed solid.
    let volume =
        ogeom::algo::volume_properties(&model, &built.shape, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    assert!(volume > 0.0 && volume.is_finite());
}

fn vertex_near(model: &Model, shape: &Shape, near: Point) -> Shape {
    explore_unique(model, shape, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .min_by(|a, b| {
            let at = |v: &Shape| {
                let p = model.node(v).unwrap().data().as_vertex().unwrap().point;
                v.transform(model.datums()).unwrap().apply(p)
            };
            at(a)
                .distance(near)
                .partial_cmp(&at(b).distance(near))
                .unwrap()
        })
        .unwrap()
}
