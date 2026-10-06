//! Offsets of fitted sheets: a filled saddle moved along its normals,
//! measured against the original by projection.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_polygon, project_on_surface};
use ogeom::core::Tolerances;
use ogeom::geom::{Continuity, Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom::math::Point;
use ogeom::offset::{FillBoundary, make_filling_n, make_thick_sheet, offset_sheet};
use ogeom::topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

fn surface(model: &Model, face: &Shape) -> SurfaceGeometry {
    let NodeData::Face(data) = model.node(face).unwrap().data() else {
        panic!("not a face");
    };
    model.geometry().surface(data.surface).unwrap().clone()
}

fn edges(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore(model, shape, Filter::OfType(ShapeType::Edge)).unwrap()
}

/// Points along an edge's curve over its range.
fn along(model: &Model, edge: &Shape, steps: u32) -> Vec<Point> {
    let NodeData::Edge(data) = model.node(edge).unwrap().data() else {
        panic!("not an edge");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        panic!("no curve");
    };
    let curve = model.geometry().curve(*curve).unwrap();
    (0..=steps)
        .map(|k| {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / f64::from(steps);
            curve.point_at(t, T).unwrap()
        })
        .collect()
}

/// The saddle filled from four straight sides: the model, the sides and
/// the face.
fn saddle() -> (Model, Vec<FillBoundary>, Shape) {
    let mut model = Model::new();
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(10.0, 0.0, 4.0),
        Point::new(10.0, 10.0, 0.0),
        Point::new(0.0, 10.0, 4.0),
    ];
    let wire = make_polygon(&mut model, &corners, true, T).unwrap().shape;
    let sides: Vec<FillBoundary> = edges(&model, &wire)
        .into_iter()
        .map(|edge| FillBoundary {
            edge,
            support: None,
            continuity: Continuity::C0,
        })
        .collect();
    let filled = make_filling_n(&mut model, &sides, &[], 1e-4, T).unwrap();
    (model, sides, filled.built.shape)
}

/// A four-sided saddle filled from four straight sides, offset by one: a
/// fitted surface the offset cannot follow exactly, fitted to a tolerance
/// of its own and reporting it.
#[test]
fn a_filled_saddle_offsets_by_one_within_a_thousandth() {
    let (mut model, sides, face) = saddle();
    let base = surface(&model, &face);

    let moved = offset_sheet(&mut model, &face, 1.0, T).unwrap().shape;
    assert_eq!(model.kind_of(&moved).unwrap(), ShapeType::Face);
    assert!(check(&model, &moved, T).unwrap().is_valid());

    // Every point of the moved surface over the face's region stands one
    // from the original, and the face's tolerance holds what the fit
    // missed by.
    let top = surface(&model, &moved);
    let ((u0, u1), (v0, v1)) = top.domain();
    let mut worst: f64 = 0.0;
    let steps = 40;
    for i in 0..=steps {
        for j in 0..=steps {
            let u = u0 + (u1 - u0) * f64::from(i) / f64::from(steps);
            let v = v0 + (v1 - v0) * f64::from(j) / f64::from(steps);
            let p = top.point_at(u, v, T).unwrap();
            let foot = project_on_surface(&base, p, 24, T).unwrap();
            worst = worst.max((foot.distance - 1.0).abs());
        }
    }
    assert!(worst <= 1e-3, "the offset strays {worst:.3e} from one");

    let reported = model.tolerance_of(&moved).unwrap().unwrap().get();
    assert!(
        reported >= worst - 1e-6,
        "the face reports {reported:.3e} against a measured {worst:.3e}"
    );

    // Four boundary edges, each the offset of one original side: every
    // point of it stands one from the base surface, its foot there on that
    // side, and each side is met by one edge.
    let bounds = edges(&model, &moved);
    assert_eq!(bounds.len(), 4);
    let sides: Vec<(Point, Point)> = sides
        .iter()
        .map(|s| {
            let ends = along(&model, &s.edge, 1);
            (ends[0], ends[1])
        })
        .collect();
    let from_side = |p: Point, (a, b): (Point, Point)| {
        let d = b - a;
        let f = ((p - a).dot(d) / d.dot(d)).clamp(0.0, 1.0);
        p.distance(a + d * f)
    };
    let mut met = vec![0; sides.len()];
    let mut edge_worst: f64 = 0.0;
    for edge in &bounds {
        let feet: Vec<(f64, Point)> = along(&model, edge, 64)
            .into_iter()
            .map(|p| {
                let foot = project_on_surface(&base, p, 24, T).unwrap();
                (foot.distance, foot.point)
            })
            .collect();
        let (index, gap) = sides
            .iter()
            .enumerate()
            .map(|(k, side)| {
                let gap = feet
                    .iter()
                    .map(|(_, q)| from_side(*q, *side))
                    .fold(0.0_f64, f64::max);
                (k, gap)
            })
            .min_by(|x, y| x.1.total_cmp(&y.1))
            .unwrap();
        assert!(gap <= 1e-3, "an edge's feet stand {gap:.3e} off every side");
        met[index] += 1;
        for (distance, _) in &feet {
            edge_worst = edge_worst.max((distance - 1.0).abs());
        }
    }
    assert!(
        edge_worst <= 1e-3,
        "an edge strays {edge_worst:.3e} from one"
    );
    assert_eq!(met, vec![1; 4]);
}

/// The saddle thickened by one: the layers fitted as the offset fits them,
/// closed into a valid solid by four side faces.
#[test]
fn a_filled_saddle_thickens_into_a_solid() {
    let (mut model, _, face) = saddle();
    let solid = make_thick_sheet(&mut model, &face, 1.0, false, T)
        .unwrap()
        .shape;
    assert_eq!(model.kind_of(&solid).unwrap(), ShapeType::Solid);
    assert!(check(&model, &solid, T).unwrap().is_valid());
    let faces = explore(&model, &solid, Filter::OfType(ShapeType::Face)).unwrap();
    assert_eq!(faces.len(), 6);
}

/// An arc in the plane z = `z` from (0, 0) to (20, 0) through (10, `bulge`).
fn bulging_arc(model: &mut Model, z: f64, bulge: f64) -> Shape {
    use ogeom::geom::CircleCurve;
    use ogeom::math::{Circle, Direction, Frame, Vector};
    let r = (100.0 + bulge * bulge) / (2.0 * bulge);
    let centre = Point::new(10.0, bulge - r, z);
    let frame = Frame::new(
        centre,
        Direction::new(Vector::new(0.0, 0.0, -1.0), T).unwrap(),
        Direction::new(Point::new(0.0, 0.0, z) - centre, T).unwrap(),
        T,
    )
    .unwrap();
    let circle = CircleCurve::new(Circle::new(frame, r, T).unwrap());
    let sweep = 2.0 * (10.0 / r).asin();
    ogeom::algo::make_edge(model, circle.into(), (0.0, sweep), T)
        .unwrap()
        .shape
}

/// The sheet between the arcs bulging 4 at z = 0 and 8 at z = 20, ruled
/// or lofted: a rational surface whose chart turns along the arcs' middle
/// knot, smooth in space but only C0 in its parameters there.
fn between_arcs(model: &mut Model, ruled: bool) -> Shape {
    use ogeom::offset::{make_loft_surface, make_ruled};
    let a = bulging_arc(model, 0.0, 4.0);
    let b = bulging_arc(model, 20.0, 8.0);
    if ruled {
        make_ruled(model, &a, &b, T).unwrap().shape
    } else {
        make_loft_surface(model, &[a, b], false, &[], false, T)
            .unwrap()
            .shape
    }
}

/// Whether the face's surface is a rational B-spline whose `u` knots hold
/// an interior knot as often as the degree: the corner line the fit has
/// to follow.
fn has_corner_line(model: &Model, face: &Shape) -> bool {
    let SurfaceGeometry::BSpline(b) = surface(model, face) else {
        return false;
    };
    let degree = b.u_knots().degree();
    b.is_rational()
        && b.u_knots().distinct().iter().any(|&(at, m)| {
            m >= degree && at > b.u_knots().domain_start() && at < b.u_knots().domain_end()
        })
}

/// The lofted sheet between two arcs of different radii offset by one:
/// every point of the moved face stands one from the source surface,
/// measured by projection, within the share of the distance a fit may miss
/// by; and thickened by three, one valid solid of six faces.
#[test]
fn a_loft_between_arcs_of_two_radii_offsets_and_thickens() {
    let mut model = Model::new();
    let sheet = between_arcs(&mut model, false);
    let face = explore(&model, &sheet, Filter::OfType(ShapeType::Face))
        .unwrap()
        .remove(0);
    assert!(has_corner_line(&model, &face));
    let base = surface(&model, &face);
    let moved = offset_sheet(&mut model, &face, 1.0, T).unwrap().shape;
    assert!(check(&model, &moved, T).unwrap().is_valid());
    let top = surface(&model, &moved);
    let ((u0, u1), (v0, v1)) = top.domain();
    let mut worst: f64 = 0.0;
    let steps = 40;
    for i in 0..=steps {
        for j in 0..=steps {
            let u = u0 + (u1 - u0) * f64::from(i) / f64::from(steps);
            let v = v0 + (v1 - v0) * f64::from(j) / f64::from(steps);
            let p = top.point_at(u, v, T).unwrap();
            let foot = project_on_surface(&base, p, 24, T).unwrap();
            worst = worst.max((foot.distance - 1.0).abs());
        }
    }
    assert!(worst <= 1e-4, "the offset strays {worst:.3e} from one");

    let solid = thickened(&mut model, &sheet, 3.0, false);
    assert_eq!(model.kind_of(&solid).unwrap(), ShapeType::Solid);
}

/// The sheet thickened, checked valid and of six faces.
fn thickened(model: &mut Model, sheet: &Shape, thickness: f64, both: bool) -> Shape {
    let solid = make_thick_sheet(model, sheet, thickness, both, T)
        .unwrap()
        .shape;
    let diagnosis = check(model, &solid, T).unwrap();
    assert!(diagnosis.is_valid(), "{thickness} {both}: {diagnosis}");
    let faces = explore(model, &solid, Filter::OfType(ShapeType::Face)).unwrap();
    assert_eq!(faces.len(), 6, "{thickness} {both}");
    solid
}

/// The ruled and the lofted sheet and the loft through a third arc between
/// them, thickened by a half to three, to one side and to both.
#[test]
#[ignore = "heavy"]
fn sheets_between_arcs_of_two_radii_thicken_at_every_distance() {
    use ogeom::offset::make_loft_surface;
    for sheet_kind in 0..3 {
        for thickness in [0.5, 1.0, 2.0, 3.0] {
            for both in [false, true] {
                let mut model = Model::new();
                let sheet = if sheet_kind < 2 {
                    between_arcs(&mut model, sheet_kind == 1)
                } else {
                    let arcs = [
                        bulging_arc(&mut model, 0.0, 4.0),
                        bulging_arc(&mut model, 10.0, 6.0),
                        bulging_arc(&mut model, 20.0, 8.0),
                    ];
                    make_loft_surface(&mut model, &arcs, false, &[], false, T)
                        .unwrap()
                        .shape
                };
                thickened(&mut model, &sheet, thickness, both);
            }
        }
    }
}
