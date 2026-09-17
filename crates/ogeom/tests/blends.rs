//! §10's tail: what a blend achieved, measured; blends between faces that
//! share no edge; edges whose envelope has no closed form; and the corner
//! where three of them meet.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The edge of `shape` whose midpoint is nearest `near`.
fn edge_near(model: &Model, shape: &Shape, near: Point) -> Shape {
    use ogeom::geom::Curve3d as _;
    // A degenerate edge — a sphere's pole — has no curve and no midpoint.
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            model
                .node(e)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|d| d.curve3d().is_some())
        })
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

/// The corner tool at every corner of the box, the three fillets in a
/// different order at each: the construction is the same whichever way the
/// corner faces and whichever edge goes first. It was not — the tool's
/// block face meets a band exactly along the arc that bounds it, and the
/// paving read that section as outside the block face by a hair at some
/// corners and inside at others, so the band split at some corners and
/// stayed whole at the rest.
#[test]
fn round_vertex_rounds_the_corner_at_any_placement() {
    let r = 3.0;
    let expected = 784.0 + 51.75 * core::f64::consts::PI;
    let fine = ogeom::mesh::Deflection::with_chord(1e-3).unwrap();
    let mut failed: Vec<(usize, String)> = Vec::new();
    for (ci, corner) in [
        (0.0, 0.0, 0.0),
        (10.0, 0.0, 0.0),
        (0.0, 10.0, 0.0),
        (10.0, 10.0, 0.0),
        (0.0, 0.0, 10.0),
        (10.0, 0.0, 10.0),
        (0.0, 10.0, 10.0),
        (10.0, 10.0, 10.0),
    ]
    .into_iter()
    .enumerate()
    {
        let mut model = Model::new();
        let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let at = Point::new(corner.0, corner.1, corner.2);
        let vertex = vertex_near(&model, &block, at);
        // The three edges' midpoints, the order rotated by the corner.
        let mut targets = [
            Point::new(5.0, corner.1, corner.2),
            Point::new(corner.0, 5.0, corner.2),
            Point::new(corner.0, corner.1, 5.0),
        ];
        targets.rotate_left(ci % 3);
        let mut solid = block;
        for target in targets {
            let edge = edge_near(&model, &solid, target);
            solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, r, T)
                .unwrap_or_else(|e| panic!("fillet at corner {ci} near {target:?}: {e}"))
                .shape;
        }
        let outcome = ogeom::fillet::round_vertex(&mut model, &solid, &vertex, r, T)
            .map_err(|e| e.to_string())
            .and_then(|rounded| {
                if !ogeom::algo::check(&model, &rounded.shape, T)
                    .unwrap()
                    .is_valid()
                {
                    return Err("not a valid solid".to_string());
                }
                ogeom::algo::volume_properties(&model, &rounded.shape, fine, T)
                    .map(|p| p.mass)
                    .map_err(|e| e.to_string())
            });
        match outcome {
            Ok(measured) => assert!(
                (measured - expected).abs() / expected < 2e-3,
                "corner {ci}: {measured} against {expected}"
            ),
            Err(e) => failed.push((ci, e)),
        }
    }
    assert!(
        failed.is_empty(),
        "placements that did not round: {failed:?}"
    );
}

/// An oblique corner: a sheared block's origin vertex, its three edges
/// filleted one after another, then the corner tool. The block is the
/// hexahedron bounded by the host planes and the three planes through the
/// ball's centre square to the edges; the patch it leaves meets its three
/// bands and three walls tangentially, and the caps at the corner are
/// consumed while the caps at the edges' far ends stand.
#[test]
fn round_vertex_rounds_an_oblique_corner() {
    let mut model = Model::new();
    let r = 2.0;
    let (a, b, c) = (
        Vector::new(20.0, 0.0, 0.0),
        Vector::new(6.0, 20.0, 0.0),
        Vector::new(3.6, 6.0, 20.0),
    );
    let block = ogeom::algo::make_parallelepiped(&mut model, Point::ORIGIN, [a, b, c], T)
        .unwrap()
        .shape;
    let vertex = vertex_near(&model, &block, Point::ORIGIN);
    let mut solid = block;
    for edge_vector in [a, b, c] {
        let edge = edge_near(&model, &solid, Point::ORIGIN + edge_vector * 0.5);
        solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, r, T)
            .unwrap()
            .shape;
    }
    let fine = ogeom::mesh::Deflection::with_chord(2e-3).unwrap();
    let before = ogeom::algo::volume_properties(&model, &solid, fine, T)
        .unwrap()
        .mass;
    let rounded = ogeom::fillet::round_vertex(&mut model, &solid, &vertex, r, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &rounded, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    // Six walls, three bands, the three far caps, and the patch.
    let faces = explore_unique(&model, &rounded, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 13, "walls, bands, far caps and the patch");
    let patch = faces
        .iter()
        .find(|f| {
            let ogeom::topo::NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            matches!(
                model.geometry().surface(data.surface),
                Some(ogeom::geom::SurfaceGeometry::Sphere(_))
            )
        })
        .expect("the corner's spherical patch");
    let contacts = ogeom::fillet::analyse_blend(&model, &rounded, patch, 15, T).unwrap();
    assert!(!contacts.is_empty());
    for contact in &contacts {
        assert!(
            contact.gap < 1e-3 && contact.tangency_error < 5e-3,
            "the patch meets its neighbour tangentially: gap {} tangency {}",
            contact.gap,
            contact.tangency_error
        );
    }
    let after = ogeom::algo::volume_properties(&model, &rounded, fine, T)
        .unwrap()
        .mass;
    assert!(
        after < before && before - after < r * r * r,
        "the corner sheds its spike and no more: {before} -> {after}"
    );
}

/// The N-support setback at a square pyramid's apex: four planes through
/// the vertex, one ball touching all four, and the corner tool's block a
/// polyhedron of eight faces. The corner is cut first and the four edges
/// then take their flush fillets one after another, each band ending on
/// the ball's rim — the other order, four fillets and then the corner,
/// dies at the third fillet, whose predecessors crash into each other at
/// the apex. The corner's volume is measured against the closed form:
/// the block, N pyramids of height `r` over the host quads, less the
/// ball's sector, whose solid angle is the apex's angular defect.
#[test]
fn round_vertex_sets_back_a_four_edge_apex() {
    let mut model = Model::new();
    let r = 1.5;
    let base_corners = [
        Point::new(-10.0, -10.0, 0.0),
        Point::new(10.0, -10.0, 0.0),
        Point::new(10.0, 10.0, 0.0),
        Point::new(-10.0, 10.0, 0.0),
    ];
    let apex = Point::new(0.0, 0.0, 15.0);
    let base = ogeom::algo::make_polygon(&mut model, &base_corners, true, T)
        .unwrap()
        .shape;
    let tip = ogeom::algo::make_vertex(&mut model, apex).shape;
    let pyramid = ogeom::offset::make_loft(&mut model, &base, &tip, T)
        .unwrap()
        .shape;
    let vertex = vertex_near(&model, &pyramid, apex);
    let fine = ogeom::mesh::Deflection::with_chord(2e-3).unwrap();
    let volume = |model: &Model, shape: &Shape| {
        ogeom::algo::volume_properties(model, shape, fine, T)
            .unwrap()
            .mass
    };
    let before = volume(&model, &pyramid);

    let rounded = ogeom::fillet::round_vertex(&mut model, &pyramid, &vertex, r, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &rounded, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert_eq!(
        explore_unique(&model, &rounded, ShapeType::Face)
            .unwrap()
            .len(),
        10,
        "five walls, the patch and four flush ends"
    );

    // The closed form, from the pyramid's own geometry: inward normals of
    // the four lateral planes, the ball's centre on the axis a radius in
    // from each, the feet of the centre on the edges and its touch points
    // on the planes.
    let edges: [Vector; 4] =
        std::array::from_fn(|k| (base_corners[k] - apex).normalized(T).unwrap());
    let inward: [Vector; 4] = std::array::from_fn(|k| {
        let n = (base_corners[(k + 1) % 4] - base_corners[k])
            .cross(apex - base_corners[k])
            .normalized(T)
            .unwrap();
        if n.dot(Point::ORIGIN - base_corners[k]) > 0.0 {
            n
        } else {
            -n
        }
    });
    let centre = apex + Vector::new(0.0, 0.0, r / inward[0].z);
    for m in &inward {
        assert!(
            (m.dot(centre - apex) - r).abs() < 1e-9,
            "one ball touches all four"
        );
    }
    let foot = |k: usize| apex + edges[k] * (centre - apex).dot(edges[k]);
    let touch = |k: usize| centre - inward[k] * r;
    let mut block = 0.0;
    let mut defect = core::f64::consts::TAU;
    for k in 0..4 {
        let (a, b, c, d) = (apex, foot(k), touch(k), foot((k + 1) % 4));
        let area = 0.5 * ((b - a).cross(c - a).magnitude() + (c - a).cross(d - a).magnitude());
        block += r * area / 3.0;
        defect -= edges[k].dot(edges[(k + 1) % 4]).acos();
    }
    let expected = block - r * r * r * defect / 3.0;
    let after_corner = volume(&model, &rounded);
    let shed = before - after_corner;
    assert!(
        (shed - expected).abs() < expected * 5e-3,
        "the apex sheds its block less the ball's sector: {shed} vs {expected}"
    );

    // The four flush fillets after the corner, each on the edge's remaining
    // run, each shedding the same volume as the others.
    let mut solid = rounded;
    let mut shed_by_band = Vec::new();
    for (k, base_corner) in base_corners.iter().enumerate() {
        let edge = edge_near(&model, &solid, foot(k).midpoint(*base_corner));
        let was = volume(&model, &solid);
        solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, r, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom::algo::check(&model, &solid, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "after fillet {k}: {:?}",
            diagnosis.problems
        );
        shed_by_band.push(was - volume(&model, &solid));
    }
    for (k, shed) in shed_by_band.iter().enumerate() {
        assert!(
            (shed - shed_by_band[0]).abs() < shed_by_band[0] * 1e-3,
            "band {k} sheds what band 0 does: {shed} vs {}",
            shed_by_band[0]
        );
    }
    let faces = explore_unique(&model, &solid, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 10, "five walls, four bands and the patch");
    // Every blend face — the patch and the four bands — meets each of its
    // neighbours tangentially: the bands their walls and the patch, the
    // patch its four bands.
    let mut blends = 0;
    for face in &faces {
        let ogeom::topo::NodeData::Face(data) = model.node(face).unwrap().data() else {
            continue;
        };
        if matches!(
            model.geometry().surface(data.surface),
            Some(ogeom::geom::SurfaceGeometry::Plane(_))
        ) {
            continue;
        }
        blends += 1;
        let contacts = ogeom::fillet::analyse_blend(&model, &solid, face, 15, T).unwrap();
        let mut tangent = 0;
        for contact in &contacts {
            // A band's flush run-out through the base is a cut, not a join:
            // its edge lies in the base plane and meets it at an angle.
            let on_base = explore_unique(&model, &contact.edge, ShapeType::Vertex)
                .unwrap()
                .iter()
                .all(|v| {
                    model
                        .node(v)
                        .unwrap()
                        .data()
                        .as_vertex()
                        .unwrap()
                        .point
                        .z
                        .abs()
                        < 1e-6
                });
            if on_base {
                continue;
            }
            tangent += 1;
            assert!(
                contact.gap < 1e-3 && contact.tangency_error < 5e-3,
                "a blend meets its neighbour tangentially: gap {} tangency {}",
                contact.gap,
                contact.tangency_error
            );
        }
        assert!(
            tangent >= 3,
            "a band meets two walls and the patch; the patch four bands"
        );
    }
    assert_eq!(blends, 5, "the patch and four bands");
}

/// A rectangular pyramid's apex: four planes, two slopes, and no ball a
/// radius in from all four at once. The tool refuses by name rather than
/// seating a ball that touches two of the faces and cuts the other two.
#[test]
fn round_vertex_refuses_an_apex_no_ball_touches() {
    let mut model = Model::new();
    let base_corners = [
        Point::new(-10.0, -5.0, 0.0),
        Point::new(10.0, -5.0, 0.0),
        Point::new(10.0, 5.0, 0.0),
        Point::new(-10.0, 5.0, 0.0),
    ];
    let apex = Point::new(0.0, 0.0, 15.0);
    let base = ogeom::algo::make_polygon(&mut model, &base_corners, true, T)
        .unwrap()
        .shape;
    let tip = ogeom::algo::make_vertex(&mut model, apex).shape;
    let pyramid = ogeom::offset::make_loft(&mut model, &base, &tip, T)
        .unwrap()
        .shape;
    let vertex = vertex_near(&model, &pyramid, apex);
    let err = ogeom::fillet::round_vertex(&mut model, &pyramid, &vertex, 1.5, T)
        .expect_err("no ball touches all four faces");
    assert!(
        err.to_string().contains("share no tangent ball"),
        "the refusal names the missing ball: {err}"
    );
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
fn an_open_seat_runs_out_through_the_wall() {
    // The bottom crease of the grooved block is an ellipse arc that meets
    // the box wall at both ends: an open seat. The band runs on past each
    // end until the ball has left the solid, and the cut trims it against
    // the wall — material comes off, the blend rides both hosts
    // tangentially, and it ends on the wall itself, not on a cap standing
    // short of it with a sliver of sharp crease behind.
    let mut model = Model::new();
    let grooved = grooved_block(&mut model);
    let before =
        ogeom::algo::volume_properties(&model, &grooved, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    let faces_before = explore_unique(&model, &grooved, ShapeType::Face)
        .unwrap()
        .len();
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
    // One new face — the band — and no caps: both ends are the wall's.
    assert_eq!(
        explore_unique(&model, &built.shape, ShapeType::Face)
            .unwrap()
            .len(),
        faces_before + 1,
        "the band is the only face the blend adds"
    );

    // The blend face is the fitted band; its rails ride the hosts
    // tangentially and every other edge of it lies on the wall.
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
    let mut on_wall = 0;
    for c in &contacts {
        assert!(c.gap < 1e-3, "a contact stands off its edge: {}", c.gap);
        if c.tangency_error < 5e-3 {
            smooth += 1;
            continue;
        }
        let NodeData::Face(d) = model.node(&c.neighbour).unwrap().data() else {
            panic!("a neighbour is a face");
        };
        let Some(ogeom::geom::SurfaceGeometry::Plane(plane)) = model.geometry().surface(d.surface)
        else {
            panic!("a non-tangent neighbour of the band is the wall, a plane");
        };
        let wall = plane.plane();
        assert!(
            wall.normal().vector().y.abs() > 0.999
                && wall.distance_to(Point::new(0.0, 20.0, 0.0)).abs() < 1e-9,
            "the band's other edges lie on the y=20 wall"
        );
        on_wall += 1;
    }
    assert_eq!(smooth, 2, "two tangent rails: {contacts:?}");
    assert!(
        on_wall >= 2,
        "the band ends on the wall at both ends: {contacts:?}"
    );
}

/// Two straight edges of a box meeting at a corner, blended together:
/// the later seat runs on through the earlier band and the cut trims the
/// two bands against each other. Each wedge removes (1 − π/4) r² per unit
/// length; the corner cell where both wedges reach is counted once, and
/// what both remove there is the cell outside both cylinders,
/// r³ (5/3 − π/2). One edge at a time stops flush instead, and keeps a
/// cap at the corner — the state the corner tool is built for.
#[test]
fn two_blends_meeting_at_a_corner_trim_each_other() {
    let (l, r) = (20.0_f64, 2.0_f64);
    let pi = core::f64::consts::PI;
    let want =
        l * l * 10.0 - (2.0 * (1.0 - pi / 4.0) * r * r * l - r * r * r * (5.0 / 3.0 - pi / 2.0));

    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (l, l, 10.0), T)
        .unwrap()
        .shape;
    let a = edge_near(&model, &block, Point::new(10.0, 20.0, 10.0));
    let b = edge_near(&model, &block, Point::new(20.0, 10.0, 10.0));
    let met = ogeom::fillet::fillet_edges(&mut model, &block, &[a.clone(), b.clone()], r, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &met, T).unwrap().is_valid());
    assert_eq!(
        explore_unique(&model, &met, ShapeType::Face).unwrap().len(),
        8,
        "six walls and two bands, no cap between them"
    );
    let mut previous = f64::INFINITY;
    for chord in [1e-3, 1e-4] {
        let fine = ogeom::mesh::Deflection::with_chord(chord).unwrap();
        let measured = ogeom::algo::volume_properties(&model, &met, fine, T)
            .unwrap()
            .mass;
        let error = (measured - want).abs() / want;
        assert!(
            error < previous,
            "refining brings it closer: {measured} vs {want}"
        );
        assert!(
            error < chord * 2.0,
            "two meeting blends against the closed form at chord {chord}: {measured} vs {want}"
        );
        previous = error;
    }

    // One edge at a time: the second stops flush at the first band, and
    // the corner cell keeps the material the meeting would have rounded.
    let first = ogeom::fillet::fillet_edge(&mut model, &block, &a, r, T).unwrap();
    let b_again = edge_near(&model, &first.shape, Point::new(20.0, 10.0, 10.0));
    let flush = ogeom::fillet::fillet_edge(&mut model, &first.shape, &b_again, r, T)
        .unwrap()
        .shape;
    assert_eq!(
        explore_unique(&model, &flush, ShapeType::Face)
            .unwrap()
            .len(),
        9,
        "flush: the second wedge's cap stands at the first band"
    );
    let fine = ogeom::mesh::Deflection::with_chord(1e-3).unwrap();
    let flush_volume = ogeom::algo::volume_properties(&model, &flush, fine, T)
        .unwrap()
        .mass;
    assert!(
        flush_volume > want + 0.5,
        "the flush corner keeps material the meeting removes"
    );
}

/// A straight seat and a marched one meeting at two corners: the box's
/// bottom edge along the wall, and the grooved block's elliptical crease
/// that runs out through that wall across it. In either order the later
/// blend's run-out walks on under the earlier band until the ball has left
/// the material, the cut trimming the two bands against each other at both
/// corners — and the two orders land on the same solid. The wall's bottom
/// edge is two pieces either side of the scoop, on one line; only the left
/// is asked for, and only the left is blended whichever goes first.
#[test]
fn a_marched_blend_meets_a_straight_blend_at_its_corners() {
    let fine = ogeom::mesh::Deflection::with_chord(2e-3).unwrap();
    let mut volumes = Vec::new();
    for straight_first in [true, false] {
        let mut model = Model::new();
        let grooved = grooved_block(&mut model);
        let before = ogeom::algo::volume_properties(&model, &grooved, fine, T)
            .unwrap()
            .mass;
        let faces_before = explore_unique(&model, &grooved, ShapeType::Face)
            .unwrap()
            .len();
        let wall_edge = edge_near(&model, &grooved, Point::new(2.0, 20.0, 0.0));
        let crease = edge_near(&model, &grooved, Point::new(10.0, 14.84, 0.0));
        let order = if straight_first {
            [wall_edge, crease]
        } else {
            [crease, wall_edge]
        };
        let met = ogeom::fillet::fillet_edges(&mut model, &grooved, &order, 1.0, T)
            .unwrap()
            .shape;
        assert!(ogeom::algo::check(&model, &met, T).unwrap().is_valid());
        assert_eq!(
            explore_unique(&model, &met, ShapeType::Face).unwrap().len(),
            faces_before + 2,
            "two bands, no cap between them, straight first {straight_first}"
        );
        // The right-hand piece of the wall's bottom edge stays sharp: the
        // side wall at x = 20 keeps its whole rectangle.
        let side = explore_unique(&model, &met, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                let b = ogeom::algo::shape_bounds(&model, f, T).unwrap();
                b.low().is_some_and(|lo| lo.x > 19.999) && b.high().is_some_and(|hi| hi.x < 20.001)
            })
            .expect("the x = 20 wall");
        let side_area = ogeom::algo::surface_properties(&model, &side, fine, T)
            .unwrap()
            .mass;
        assert!(
            (side_area - 200.0).abs() < 1e-3,
            "the unrequested right piece is not blended: {side_area}"
        );
        let after = ogeom::algo::volume_properties(&model, &met, fine, T)
            .unwrap()
            .mass;
        assert!(after < before && after > before * 0.9);
        volumes.push(after);
    }
    assert!(
        (volumes[0] - volumes[1]).abs() < 1e-2,
        "both orders round the same material: {volumes:?}"
    );
}

#[test]
fn two_seam_split_blends_meet_cap_to_cap_in_either_order() {
    // The top crease is split by the drum's seam into two arcs sharing a
    // vertex mid-scoop. Each rounds as a capped blend ending in the arc's
    // own section plane at the seam vertex; the second blend's cap meets
    // the first's in that plane, the two bands meet along the shared arc,
    // and both caps are consumed. Whichever arc goes first, the result is
    // the same closed solid with two bands and no cap.
    let mut volumes = Vec::new();
    for order in [
        [Point::new(7.3, 7.7, 10.0), Point::new(12.7, 7.7, 10.0)],
        [Point::new(12.7, 7.7, 10.0), Point::new(7.3, 7.7, 10.0)],
    ] {
        let mut model = Model::new();
        let grooved = grooved_block(&mut model);
        let faces_before = explore_unique(&model, &grooved, ShapeType::Face)
            .unwrap()
            .len();
        let before =
            ogeom::algo::volume_properties(&model, &grooved, ogeom::mesh::Deflection::default(), T)
                .unwrap()
                .mass;
        let first_arc = edge_near(&model, &grooved, order[0]);
        let first = ogeom::fillet::fillet_edge(&mut model, &grooved, &first_arc, 1.0, T).unwrap();
        let second_arc = edge_near(&model, &first.shape, order[1]);
        let second =
            ogeom::fillet::fillet_edge(&mut model, &first.shape, &second_arc, 1.0, T).unwrap();
        for shell in explore_unique(&model, &second.shape, ShapeType::Shell).unwrap() {
            assert!(ogeom::algo::is_shell_closed(&model, &shell).unwrap());
        }
        assert_eq!(
            explore_unique(&model, &second.shape, ShapeType::Face)
                .unwrap()
                .len(),
            faces_before + 2,
            "two bands, no caps"
        );
        let after = ogeom::algo::volume_properties(
            &model,
            &second.shape,
            ogeom::mesh::Deflection::default(),
            T,
        )
        .unwrap()
        .mass;
        assert!(after < before && after > before * 0.9);
        volumes.push(after);
    }
    assert!(
        (volumes[0] - volumes[1]).abs() < 1e-2,
        "the order does not change the solid: {volumes:?}"
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

#[test]
fn two_disjoint_run_out_blends_coexist_on_one_solid() {
    // The first capped blend lands on the bottom crease; the second on the
    // top crease's far seam-half, nowhere near the first. Sequential
    // marched blends must not disturb each other's wounds.
    let mut model = Model::new();
    let grooved = grooved_block(&mut model);
    let first_arc = edge_near(&model, &grooved, Point::new(10.0, 14.84, 0.0));
    let first = ogeom::fillet::fillet_edge(&mut model, &grooved, &first_arc, 1.0, T).unwrap();
    let second_arc = edge_near(&model, &first.shape, Point::new(7.3, 7.7, 10.0));
    let second = ogeom::fillet::fillet_edge(&mut model, &first.shape, &second_arc, 1.0, T).unwrap();
    let volume = ogeom::algo::volume_properties(
        &model,
        &second.shape,
        ogeom::mesh::Deflection::default(),
        T,
    )
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
