//! A face whose boundary, read edge by edge in its chart, jumps a whole
//! period at one joint and back at another: an edge on the closure of a
//! closed chart stated on the other branch from the edges it meets. The
//! ring is one closed curve on the surface, and a boolean through the face
//! reads it as one.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_box, make_cone, make_cylinder, make_sphere, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::{Curve2d as _, Curve3d as _, PlanarCurve};
use ogeom::math::{Direction, Frame, Point, Transform2, Vector, Vector2};
use ogeom::mesh::Deflection;
use ogeom::topo::{EdgeRepr, Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const TAU: f64 = core::f64::consts::TAU;

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

fn block(model: &mut Model, corner: (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
    let at = Frame::new(
        Point::new(corner.0, corner.1, corner.2),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    make_box(model, at, size, T).unwrap().shape
}

/// A bar 4 long along x on a 2 by 2 square, with a quarter round of radius
/// 1 notched out along its edge at y = z = 2: the cylinder's chart starts
/// at `-y` and turns towards `-z` or `+z` as `axis` points along `+x` or
/// `-x`, so the wall spans the chart's first quarter or its last, and one
/// of its straight edges lies on the chart's closure either way. That
/// edge's pcurve on the wall is then moved a whole period to the other
/// line of the closure.
fn folded_notched_bar(model: &mut Model, axis: Direction) -> Shape {
    let bar = block(model, (0.0, 0.0, 0.0), (4.0, 2.0, 2.0));
    let from = if axis == Direction::X { -1.0 } else { 5.0 };
    let pin = make_cylinder(
        model,
        Frame::new(
            Point::new(from, 2.0, 2.0),
            axis,
            Direction::new(ogeom::math::Vector::new(0.0, -1.0, 0.0), T).unwrap(),
            T,
        )
        .unwrap(),
        1.0,
        6.0,
        T,
    )
    .unwrap()
    .shape;
    let notched = ogeom::boolean::cut(model, &bar, &pin, T).unwrap().shape;
    let mut folded = 0;
    for face in explore_unique(model, &notched, ShapeType::Face).unwrap() {
        let surface = model.node(&face).unwrap().data().as_face().unwrap().surface;
        if !matches!(
            model.geometry().surface(surface),
            Some(ogeom::geom::SurfaceGeometry::Cylinder(_))
        ) {
            continue;
        }
        // Each wall edge's column in the chart as the face reads it, and
        // where the wall's arcs reach: the closure's line beside them is
        // the one the ring walks.
        let mut columns = Vec::new();
        let mut arcs_at_zero = false;
        for edge in explore_unique(model, &face, ShapeType::Edge).unwrap() {
            let data = model.node(&edge).unwrap().data().as_edge().unwrap();
            let (curve, range, location) = match data.pcurve_for(surface, edge.location()) {
                Some(EdgeRepr::PCurve {
                    curve,
                    range,
                    location,
                    ..
                })
                | Some(EdgeRepr::Seam {
                    forward: curve,
                    range,
                    location,
                    ..
                }) => (*curve, *range, location.clone()),
                _ => continue,
            };
            let stated = model.geometry().pcurve(curve).unwrap().clone();
            let (a, b) = (
                stated.point_at(range.0, T).unwrap(),
                stated.point_at(range.1, T).unwrap(),
            );
            if (a.x - b.x).abs() > 1e-9 {
                arcs_at_zero |= a.x.min(b.x).abs() < 1e-9;
                continue;
            }
            let on_closure = a.x.abs() < 1e-9 || (a.x - TAU).abs() < 1e-9;
            if on_closure {
                columns.push((edge, stated, range, location, a.x));
            }
        }
        for (edge, stated, range, location, u) in columns {
            let far = if arcs_at_zero { TAU } else { 0.0 };
            let moved: PlanarCurve = stated
                .transformed(&Transform2::translation(Vector2::new(far - u, 0.0)), T)
                .unwrap();
            let moved = model.geometry_mut().add_pcurve(moved);
            if let Some(node) = model.node_mut(&edge)
                && let NodeData::Edge(data) = node.data_mut()
            {
                data.representations
                    .retain(|r| r.surface() != Some(surface));
                data.representations.push(EdgeRepr::PCurve {
                    curve: moved,
                    range,
                    surface,
                    location,
                });
            }
            folded += 1;
        }
    }
    assert_eq!(
        folded, 1,
        "one straight edge of the wall lies on the closure"
    );
    notched
}

#[test]
fn a_bar_with_a_folded_wall_ring_cuts_through_the_wall() {
    // The bar is the 2 by 2 square less a quarter disc of radius 1, over
    // its length: a slab 1 thick across it, through the wall, leaves 3 of
    // its 4.
    let profile = 4.0 - core::f64::consts::FRAC_PI_4;
    for axis in [Direction::X, -Direction::X] {
        let mut model = Model::with_tolerances(T);
        let bar = folded_notched_bar(&mut model, axis);
        let found = check(&model, &bar, T).unwrap();
        assert!(found.is_valid(), "{axis:?}: {found}");
        let slab = block(&mut model, (1.5, -1.0, -1.0), (1.0, 4.0, 4.0));
        let cut = ogeom::boolean::cut(&mut model, &bar, &slab, T)
            .unwrap_or_else(|e| panic!("{axis:?}: {e}"))
            .shape;
        let found = check(&model, &cut, T).unwrap();
        assert!(found.is_valid(), "{axis:?}: {found}");
        let got = volume(&model, &cut);
        assert!(
            (got - 3.0 * profile).abs() <= 1e-6,
            "{axis:?}: {got} against {}",
            3.0 * profile
        );
    }
}

/// A solid with its closed charts' seams turned to 225 degrees, cut down
/// to the quadrant x > 0, y > 0, where no seam reaches: every face on a
/// closed chart is seam-free, its rings stated on the chart's own branch.
fn quadrant_of(model: &mut Model, kind: &str) -> Shape {
    let along = |x: f64, y: f64, z: f64| Direction::new(Vector::new(x, y, z), T).unwrap();
    let frame = |o: Point, z: Direction, x: Direction| Frame::new(o, z, x, T).unwrap();
    let seam = along(-1.0, -1.0, 0.0);
    let quadrant = block(model, (0.0, 0.0, -10.0), (10.0, 10.0, 20.0));
    let whole = match kind {
        "cone" => make_cone(
            model,
            frame(Point::ORIGIN, Direction::Z, seam),
            3.0,
            0.0,
            4.0,
            T,
        ),
        "sphere" => make_sphere(model, frame(Point::ORIGIN, Direction::Z, seam), 3.0, T),
        _ => make_cylinder(model, frame(Point::ORIGIN, Direction::Z, seam), 3.0, 6.0, T),
    }
    .unwrap()
    .shape;
    let part = ogeom::boolean::common(model, &whole, &quadrant, T)
        .unwrap()
        .shape;
    if kind != "drilled wall" {
        return part;
    }
    let a = 0.6_f64;
    let pin = make_cylinder(
        model,
        frame(
            Point::new(0.0, 0.0, 3.0),
            along(a.cos(), a.sin(), 0.0),
            Direction::Z,
        ),
        0.5,
        10.0,
        T,
    )
    .unwrap()
    .shape;
    ogeom::boolean::cut(model, &part, &pin, T).unwrap().shape
}

/// The plain pcurves of the first face of `solid` (its face on a closed
/// chart), each a case: the face's index of the edge and the period.
fn restatable(model: &Model, solid: &Shape) -> Vec<(usize, f64)> {
    let face = explore_unique(model, solid, ShapeType::Face).unwrap()[0].clone();
    let surface = model.node(&face).unwrap().data().as_face().unwrap().surface;
    let mut cases = Vec::new();
    for (k, edge) in explore_unique(model, &face, ShapeType::Edge)
        .unwrap()
        .iter()
        .enumerate()
    {
        let data = model.node(edge).unwrap().data().as_edge().unwrap();
        if data.curve3d().is_some()
            && matches!(
                data.pcurve_for(surface, edge.location()),
                Some(EdgeRepr::PCurve { .. })
            )
        {
            cases.extend([(k, TAU), (k, -TAU)]);
        }
    }
    cases
}

/// The `k`-th edge of the first face of `solid`, its pcurve on that face
/// moved `by` along `u`.
fn restate(model: &mut Model, solid: &Shape, k: usize, by: f64) {
    let face = explore_unique(model, solid, ShapeType::Face).unwrap()[0].clone();
    let surface = model.node(&face).unwrap().data().as_face().unwrap().surface;
    let edge = explore_unique(model, &face, ShapeType::Edge).unwrap()[k].clone();
    let Some(EdgeRepr::PCurve {
        curve,
        range,
        location,
        ..
    }) = model
        .node(&edge)
        .unwrap()
        .data()
        .as_edge()
        .unwrap()
        .pcurve_for(surface, edge.location())
        .cloned()
    else {
        unreachable!("a plain pcurve")
    };
    let moved: PlanarCurve = model
        .geometry()
        .pcurve(curve)
        .unwrap()
        .transformed(&Transform2::translation(Vector2::new(by, 0.0)), T)
        .unwrap();
    let moved = model.geometry_mut().add_pcurve(moved);
    if let Some(node) = model.node_mut(&edge)
        && let NodeData::Edge(data) = node.data_mut()
    {
        data.representations
            .retain(|r| r.surface() != Some(surface));
        data.representations.push(EdgeRepr::PCurve {
            curve: moved,
            range,
            surface,
            location,
        });
    }
}

/// The cut and the common of `solid` with the half-spaces x > 1.137 and
/// y > 1.113: each one's volume.
fn halves(model: &mut Model, solid: &Shape) -> Vec<f64> {
    let mut out = Vec::new();
    for corner in [(1.137, -50.0, -50.0), (-50.0, 1.113, -50.0)] {
        let tool = block(model, corner, (100.0, 100.0, 100.0));
        for common in [false, true] {
            let made = if common {
                ogeom::boolean::common(model, solid, &tool, T)
            } else {
                ogeom::boolean::cut(model, solid, &tool, T)
            }
            .unwrap()
            .shape;
            let found = check(model, &made, T).unwrap();
            assert!(found.is_valid(), "{found}");
            out.push(volume(model, &made));
        }
    }
    out
}

#[test]
fn a_restated_edge_cuts_as_the_solid_stated_on_one_branch_does() {
    // Each plain pcurve of the face on the closed chart, moved a period
    // either way, leaves the same curve on the surface: the cuts and
    // commons through it measure as the unmoved solid's do. A cone's
    // sector carries an apex edge whose chart line anchors its ring, a
    // sphere's lune two pole edges, and a wall a hole drilled through it.
    for kind in ["cone", "sphere", "wall", "drilled wall"] {
        let mut model = Model::with_tolerances(T);
        let solid = quadrant_of(&mut model, kind);
        let want = halves(&mut model, &solid);
        let cases = restatable(&model, &solid);
        assert!(!cases.is_empty(), "{kind}");
        for (k, by) in cases {
            let mut model = Model::with_tolerances(T);
            let solid = quadrant_of(&mut model, kind);
            restate(&mut model, &solid, k, by);
            assert!(
                check(&model, &solid, T).unwrap().is_valid(),
                "{kind} {k} {by}"
            );
            let got = halves(&mut model, &solid);
            for (g, w) in got.iter().zip(&want) {
                assert!(
                    (g - w).abs() <= 1e-9 * w.max(1.0),
                    "{kind}: edge {k} moved {by}: {got:?} against {want:?}"
                );
            }
        }
    }
}

/// The drilled wall's face on the cylinder, alone as a sheet, with the two
/// edges its hole is made of joined into one closed edge: the section the
/// drill cut over its whole range, on the pcurve the two pieces were cut
/// from, closing at the vertex where that range starts.
fn wall_with_a_single_edge_hole(model: &mut Model) -> Shape {
    let solid = quadrant_of(model, "drilled wall");
    let face = explore_unique(model, &solid, ShapeType::Face).unwrap()[0].clone();
    let surface = model.node(&face).unwrap().data().as_face().unwrap().surface;
    let hole = model
        .ordered_children_of(&model.children_of(&face).unwrap()[1])
        .unwrap();
    assert_eq!(hole.len(), 2);
    let data: Vec<_> = hole
        .iter()
        .map(|e| model.node(e).unwrap().data().as_edge().unwrap().clone())
        .collect();
    let piece = |d: &ogeom::topo::EdgeData| {
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = d.curve3d().cloned() else {
            unreachable!("a curve")
        };
        let Some(EdgeRepr::PCurve { curve: pcurve, .. }) = d
            .pcurve_for(surface, &ogeom::topo::Location::identity())
            .cloned()
        else {
            unreachable!("a pcurve on the wall")
        };
        (
            model.geometry().curve(curve).unwrap().clone(),
            model.geometry().pcurve(pcurve).unwrap().clone(),
            range,
        )
    };
    // The pieces stand on copies of one section and one pcurve over the
    // whole range: each reads as the other.
    let (first, second) = (piece(&data[0]), piece(&data[1]));
    for t in [0.1, 0.5, 0.9] {
        assert!(
            first
                .0
                .point_at(t, T)
                .unwrap()
                .distance(second.0.point_at(t, T).unwrap())
                < 1e-12
        );
        assert!(
            first
                .1
                .point_at(t, T)
                .unwrap()
                .distance(second.1.point_at(t, T).unwrap())
                < 1e-12
        );
    }
    let starts = if first.2.0 == 0.0 { 0 } else { 1 };
    let start =
        model.children_of(&hole[starts]).unwrap()[0].oriented(ogeom::topo::Orientation::Forward);
    let range = (0.0, 1.0);
    let (curve, pcurve) = (first.0, first.1);
    let edge = ogeom::algo::make_edge_between(model, curve, range, &start, &start, T)
        .unwrap()
        .shape;
    let pcurve = model.geometry_mut().add_pcurve(pcurve);
    let tolerance = data.iter().map(|d| d.tolerance.get()).fold(0.0, f64::max);
    model
        .widen(&edge, ogeom::core::Tolerance::new(tolerance).unwrap())
        .unwrap();
    if let Some(node) = model.node_mut(&edge)
        && let NodeData::Edge(data) = node.data_mut()
    {
        data.representations.push(EdgeRepr::PCurve {
            curve: pcurve,
            range,
            surface,
            location: ogeom::topo::Location::identity(),
        });
    }
    let mut reshape = ogeom::heal::Reshape::new();
    reshape.replace(
        &hole[0],
        if hole[0].orientation() == ogeom::topo::Orientation::Reversed {
            edge.reversed()
        } else {
            edge
        },
    );
    reshape.remove(&hole[1]);
    let face = reshape.apply(model, &face).unwrap().shape;
    ogeom::algo::build::make_shell(model, &[face])
        .unwrap()
        .shape
}

/// The area of the cut and the common of `sheet` with the half-spaces
/// x > 1.137 and y > 1.113.
fn sheet_halves(model: &mut Model, sheet: &Shape) -> Vec<f64> {
    let mut out = Vec::new();
    for corner in [(1.137, -50.0, -50.0), (-50.0, 1.113, -50.0)] {
        let tool = block(model, corner, (100.0, 100.0, 100.0));
        for common in [false, true] {
            let made = if common {
                ogeom::boolean::common(model, sheet, &tool, T)
            } else {
                ogeom::boolean::cut(model, sheet, &tool, T)
            }
            .unwrap()
            .shape;
            let found = check(model, &made, T).unwrap();
            assert!(found.is_usable(), "{found}");
            out.push(
                ogeom::algo::surface_properties(
                    model,
                    &made,
                    Deflection::with_chord(1e-3).unwrap(),
                    T,
                )
                .unwrap()
                .mass,
            );
        }
    }
    out
}

#[test]
fn a_single_edge_hole_is_read_on_its_outer_ring_s_branch() {
    // A quarter wall's sheet with a hole of one closed edge, its pcurve
    // moved a period either way: the curve on the surface is the same, so
    // the sheet cut and kept by two half-spaces through the hole measures
    // as the unmoved sheet does.
    let mut model = Model::with_tolerances(T);
    let sheet = wall_with_a_single_edge_hole(&mut model);
    assert!(check(&model, &sheet, T).unwrap().is_usable());
    let want = sheet_halves(&mut model, &sheet);
    for by in [TAU, -TAU] {
        let mut model = Model::with_tolerances(T);
        let sheet = wall_with_a_single_edge_hole(&mut model);
        let face = explore_unique(&model, &sheet, ShapeType::Face).unwrap()[0].clone();
        let hole = model
            .ordered_children_of(&model.children_of(&face).unwrap()[1])
            .unwrap();
        assert_eq!(hole.len(), 1);
        let k = explore_unique(&model, &face, ShapeType::Edge)
            .unwrap()
            .iter()
            .position(|e| e.is_same(&hole[0]))
            .unwrap();
        restate(&mut model, &sheet, k, by);
        assert!(check(&model, &sheet, T).unwrap().is_usable(), "{by}");
        let got = sheet_halves(&mut model, &sheet);
        for (g, w) in got.iter().zip(&want) {
            assert!(
                (g - w).abs() <= 1e-9 * w.max(1.0),
                "hole moved {by}: {got:?} against {want:?}"
            );
        }
    }
}
