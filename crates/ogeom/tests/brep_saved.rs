//! Shapes as an application saves them in the `.brep` text format: written
//! straight from the shape with no drawing header, placements kept as
//! chains of elementary ones, charts bound to the placement an edge has in
//! its face, wires listing their edges as a set, surfaces on left-handed
//! seats. Each fixture under `saved_shapes/` is a small shape written for
//! these tests, and every expectation is the closed form of the shape it
//! describes.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::{Curve3d as _, SurfaceGeometry};
use ogeom::math::Point;
use ogeom::mesh::Deflection;
use ogeom::topo::{EdgeRepr, Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

const PADDED_DISC: &str = include_str!("saved_shapes/padded_disc.brp");
const PADDED_ELL: &str = include_str!("saved_shapes/padded_ell.brp");
const PLACED_BLOCK: &str = include_str!("saved_shapes/placed_block.brp");
const MIRRORED_DRILL: &str = include_str!("saved_shapes/mirrored_drill.brp");
const BORED_BLOCK: &str = include_str!("saved_shapes/bored_block.brp");

/// The text as the drawing tool would write it: the content-type line
/// first, which the saved shapes leave out.
fn with_header(text: &str) -> String {
    format!("DBRep_DrawableShape\n{text}")
}

fn volume(model: &Model, shape: &Shape) -> ogeom::algo::MassProperties {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T).unwrap()
}

fn near(a: Point, b: Point) -> bool {
    a.distance(b) < 1e-9
}

/// The middle of an edge's curve, where the edge is placed.
fn middle(model: &Model, edge: &Shape) -> Point {
    let data = model.node(edge).unwrap().data().as_edge().unwrap();
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        panic!("an edge with no curve");
    };
    let at = model
        .geometry()
        .curve(*curve)
        .unwrap()
        .point_at(f64::midpoint(range.0, range.1), T)
        .unwrap();
    edge.transform(model.datums()).unwrap().apply(at)
}

#[test]
fn a_shape_saved_without_the_drawing_header_reads() {
    // The block written straight from the shape opens with a blank line and
    // the version line; it is the same file as one with the header.
    assert!(!PLACED_BLOCK.starts_with("DBRep_DrawableShape"));
    let (model, shape) = ogeom::io::brep::read(PLACED_BLOCK, T).unwrap();
    // A 4 by 2 by 1 block. Planar faces integrate exactly, so only
    // rounding stands between the sum and eight.
    let mass = volume(&model, &shape).mass;
    assert!((mass - 8.0).abs() < 1e-9, "a 4 x 2 x 1 block: {mass}");
}

#[test]
fn a_seam_whose_continuity_runs_into_its_number_reads() {
    // The disc's seam is written `1 2CN`: chart 1, chart 2, continuity CN,
    // the last two with no space between them.
    assert!(PADDED_DISC.contains(" 2CN "));
    let (model, shape) = ogeom::io::brep::read(&with_header(PADDED_DISC), T).unwrap();
    let seams = explore_unique(&model, &shape, ShapeType::Edge)
        .unwrap()
        .iter()
        .filter(|edge| {
            let data = model.node(edge).unwrap().data().as_edge().unwrap();
            data.representations
                .iter()
                .any(|r| matches!(r, EdgeRepr::Seam { .. }))
        })
        .count();
    assert_eq!(seams, 1, "the cylinder's one seam, with both its charts");
}

#[test]
fn a_composed_placement_applies_its_first_factor_first() {
    // The block, corner at the origin, is placed by a chain of two: a
    // quarter turn about x, then a shift of (3, -2, 7). Its middle (2, 1,
    // 0.5) turns to (2, -0.5, 1) and shifts to (5, -2.5, 8); the shift
    // applied first would put it at (5, -7.5, -1).
    let (model, shape) = ogeom::io::brep::read(&with_header(PLACED_BLOCK), T).unwrap();
    let centre = volume(&model, &shape).centre;
    assert!(
        near(centre, Point::new(5.0, -2.5, 8.0)),
        "turned, then shifted: {centre:?}"
    );
}

#[test]
fn a_chart_bound_to_a_placed_edge_bounds_its_face_there() {
    // The cylinder's two rims are one circle placed twice, and each rim's
    // chart on the cylinder is bound to its own placement: read without
    // them, both rims would take the first chart and the wall would close
    // on one rim. The continuity is spaced here, and the wall measured
    // alone, to read the placements alone.
    let text = with_header(PADDED_DISC).replace(" 2CN ", " 2 CN ");
    let (model, shape) = ogeom::io::brep::read(&text, T).unwrap();
    let wall = explore_unique(&model, &shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|face| {
            let data = model.node(face).unwrap().data().as_face().unwrap();
            matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::Cylinder(_))
            )
        })
        .unwrap();
    let area = ogeom::algo::surface_properties(&model, &wall, Deflection::default(), T)
        .unwrap()
        .mass;
    let exact = 2.0 * core::f64::consts::PI * 10.0 * 5.0;
    // The wall between the rims is a chart rectangle, integrated in closed
    // form: a millionth is rounding's share many times over.
    assert!(
        (area - exact).abs() < exact * 1e-6,
        "a wall of radius 10 between rims 5 apart: {area} against {exact}"
    );
}

#[test]
fn a_saved_cylinder_reads_as_itself_numbered_as_the_file_numbers_it() {
    let (model, shape, numbering) = ogeom::io::brep::read_numbered(PADDED_DISC, T).unwrap();
    let mass = volume(&model, &shape).mass;
    let exact = core::f64::consts::PI * 10.0 * 10.0 * 5.0;
    // The wall is a cylinder and the caps planes; the closed forms carry
    // them to a few ulps of the sum, far inside a millionth.
    assert!(
        (mass - exact).abs() < exact * 1e-6,
        "a 10 mm radius, 5 mm long cylinder: {mass} against {exact}"
    );

    // The file numbers the wall first, then the cap at the placement
    // (turn, then shift by (3, -2, 7)) and the cap moved 5 further down y.
    let faces = explore_unique(&model, &shape, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 3, "a wall and two caps");
    let surface = |face: &Shape| {
        let data = model.node(face).unwrap().data().as_face().unwrap();
        model.geometry().surface(data.surface).unwrap().clone()
    };
    assert!(matches!(
        surface(&numbering.faces[0]),
        SurfaceGeometry::Cylinder(_)
    ));
    let centres: Vec<Point> = numbering.faces[1..]
        .iter()
        .map(|face| {
            ogeom::algo::surface_properties(&model, face, Deflection::default(), T)
                .unwrap()
                .centre
        })
        .collect();
    assert!(near(centres[0], Point::new(3.0, -2.0, 7.0)), "{centres:?}");
    assert!(near(centres[1], Point::new(3.0, -7.0, 7.0)), "{centres:?}");
    let key = |s: &Shape| (s.node(), s.location().clone());
    assert_eq!(
        numbering.faces.iter().map(key).collect::<Vec<_>>(),
        faces.iter().map(key).collect::<Vec<_>>(),
        "the shape's faces come in the file's order"
    );

    // And its edges: the seam, then the rim at the placement, then the rim
    // further down, as the wall's wire lists them. The wall's wire is
    // walked end to end in the shape, which meets the rims the other way
    // round, so it is the numbering that keeps the file's order.
    let middles: Vec<Point> = numbering.edges.iter().map(|e| middle(&model, e)).collect();
    assert_eq!(middles.len(), 3);
    assert!(near(middles[0], Point::new(13.0, -4.5, 7.0)), "{middles:?}");
    assert!(near(middles[1], Point::new(-7.0, -2.0, 7.0)), "{middles:?}");
    assert!(near(middles[2], Point::new(-7.0, -7.0, 7.0)), "{middles:?}");
    let mut edges: Vec<_> = explore_unique(&model, &shape, ShapeType::Edge)
        .unwrap()
        .iter()
        .map(key)
        .collect();
    let mut numbered: Vec<_> = numbering.edges.iter().map(key).collect();
    edges.sort_by_key(|k| format!("{k:?}"));
    numbered.sort_by_key(|k| format!("{k:?}"));
    assert_eq!(numbered, edges, "the numbering names the shape's own edges");
}

#[test]
fn saved_shapes_measure_their_closed_forms() {
    let pi = core::f64::consts::PI;
    // The drill's point is 118 degrees across, so its cone rises
    // r / tan(59 degrees) over the bore's 2.5 mm radius.
    let (r, h) = (2.5_f64, 12.0);
    let tip = r / 59.0_f64.to_radians().tan();
    let cases = [
        // An L of area 150, padded 7: side walls that list their edges out
        // of walking order, and caps that are one face placed twice.
        ("padded ell", PADDED_ELL, 150.0 * 7.0),
        // Mirrored, so its wall and its point sit on left-handed seats,
        // the point's chart running along its ruling.
        (
            "mirrored drill",
            MIRRORED_DRILL,
            pi * r * r * (h + tip / 3.0),
        ),
        // A 30 by 20 by 10 block less a bore of radius 4 right through it.
        (
            "bored block",
            BORED_BLOCK,
            30.0 * 20.0 * 10.0 - pi * 4.0 * 4.0 * 10.0,
        ),
    ];
    for (what, text, exact) in cases {
        let (model, shape) = ogeom::io::brep::read(text, T).unwrap();
        let mass = volume(&model, &shape).mass;
        // Planes, cylinders and cones all integrate in closed form; a
        // millionth is rounding's share many times over.
        assert!(
            (mass - exact).abs() < exact * 1e-6,
            "{what}: {mass} against {exact}"
        );
        let check = ogeom::algo::check(&model, &shape, T).unwrap();
        assert!(check.is_valid(), "{what}: {check:?}");
    }
}
