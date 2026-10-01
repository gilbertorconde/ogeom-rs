//! Filleting a tangent chain of edges in one call.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::CircleCurve;
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Filter, ShapeType, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &ogeom_topo::Model, shape: &ogeom_topo::Shape, chord: f64) -> f64 {
    ogeom_algo::volume_properties(
        model,
        shape,
        ogeom_mesh::Deflection {
            chord,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

/// A stadium prism: a rectangle with semicircular ends, extruded up.
fn stadium(model: &mut ogeom_topo::Model, length: f64, r: f64, height: f64) -> ogeom_topo::Shape {
    // Four corner vertices, shared between neighbouring edges: a wire chains
    // through vertex objects, not coincident coordinates.
    let tl = ogeom_algo::make_vertex(model, Point::new(0.0, r, 0.0)).shape;
    let tr = ogeom_algo::make_vertex(model, Point::new(length, r, 0.0)).shape;
    let br = ogeom_algo::make_vertex(model, Point::new(length, -r, 0.0)).shape;
    let bl = ogeom_algo::make_vertex(model, Point::new(0.0, -r, 0.0)).shape;

    // Each arc runs (0, pi) on a frame whose x points at its own start, so
    // every window is canonical.
    let arc = |model: &mut ogeom_topo::Model,
               centre: Point,
               x: Direction,
               from: &ogeom_topo::Shape,
               to: &ogeom_topo::Shape|
     -> ogeom_topo::Shape {
        let frame = Frame::new(centre, Direction::Z, x, T).unwrap();
        let circle = Circle::new(frame, r, T).unwrap();
        let curve = ogeom_geom::Curve::Circle(CircleCurve::new(circle));
        ogeom_algo::make_edge_between(model, curve, (0.0, core::f64::consts::PI), from, to, T)
            .unwrap()
            .shape
    };
    let seg = |model: &mut ogeom_topo::Model,
               from: (&ogeom_topo::Shape, Point),
               to: (&ogeom_topo::Shape, Point)|
     -> ogeom_topo::Shape {
        let line = ogeom_geom::LineCurve::segment(from.1, to.1, T).unwrap();
        let curve = ogeom_geom::Curve::Line(line);
        let domain = ogeom_geom::Curve3d::domain(&curve);
        ogeom_algo::make_edge_between(model, curve, domain, from.0, to.0, T)
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
        Direction::new(ogeom_math::Vector::new(0.0, -1.0, 0.0), T).unwrap(),
        &br,
        &tr,
    );
    let bottom = seg(
        model,
        (&br, Point::new(length, -r, 0.0)),
        (&bl, Point::new(0.0, -r, 0.0)),
    );
    let left = arc(model, Point::new(0.0, 0.0, 0.0), Direction::Y, &tl, &bl);
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face_with_pcurves(
        model,
        plane.into(),
        &[vec![left, bottom.reversed(), right, top.reversed()]],
        T,
    )
    .unwrap()
    .shape;
    ogeom_algo::make_prism(model, &face, Vector::new(0.0, 0.0, height), T)
        .unwrap()
        .shape
}

#[test]
fn a_tangent_chain_of_edges_fillets_in_one_call() {
    // The stadium's top rim: two straight runs and two semicircular ends,
    // tangent all the way round. One call blends the four; at each junction
    // the neighbouring wedges' end caps stand in one plane with one
    // cross-section, and the melt joins the blends without a seam.
    let mut model = ogeom_topo::Model::new();
    let (length, r, height, blend) = (10.0, 5.0, 4.0, 1.0);
    let solid = stadium(&mut model, length, r, height);
    let before = volume(&model, &solid, 1e-3);

    // The rim: the four top edges, ordered around the loop.
    let rim: Vec<ogeom_topo::Shape> = explore_unique(&model, &solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            explore(&model, e, Filter::OfType(ShapeType::Vertex))
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
        .collect();
    assert_eq!(rim.len(), 4, "the stadium's top rim has four edges");

    let result = ogeom_fillet::fillet_edges(&mut model, &solid, &rim, blend, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    for edge in &rim {
        assert!(result.history.is_deleted(edge), "every rim edge is gone");
    }

    // The removed ring: the fillet's cross-section swept along the rim. On
    // the straight runs that is exactly (1 - pi/4) r^2 per unit length; on
    // the arcs, Pappus moves it by the cross-section's centroid, whose
    // offset from the rim is known in closed form.
    let pi = core::f64::consts::PI;
    let section = (1.0 - pi / 4.0) * blend * blend;
    let straight = 2.0 * length * section;
    // Centroid of the region between a square and its inscribed quarter
    // disc, measured from the rim corner, resolved radially inward.
    let centroid = (10.0 - 3.0 * pi) / (12.0 - 3.0 * pi) * blend;
    let arcs = 2.0 * pi * (r - centroid) * section;
    let expected = before - straight - arcs;
    let measured = volume(&model, &result.shape, 1e-3);
    assert!(
        (measured - expected).abs() < 0.05,
        "rounded stadium volume {measured} against {expected}"
    );

    // The blends themselves: cylinders along the straights, tori round the
    // ends, and no leftover cap faces at the tangent junctions.
    let mut cylinders = 0;
    let mut tori = 0;
    for f in explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap() {
        match model
            .node(&f)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
        {
            Some(ogeom_geom::SurfaceGeometry::Cylinder(_)) => cylinders += 1,
            Some(ogeom_geom::SurfaceGeometry::Torus(_)) => tori += 1,
            _ => {}
        }
    }
    assert!(
        cylinders >= 4,
        "wall drums and straight blends: {cylinders}"
    );
    assert_eq!(tori, 2, "one torus blend per rounded end: {tori}");
}

/// The edge of `shape` whose middle stands at `at`.
fn edge_through(
    model: &ogeom_topo::Model,
    shape: &ogeom_topo::Shape,
    at: Point,
) -> ogeom_topo::Shape {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (a, b) = ogeom_algo::edge_vertices(model, e).unwrap().unwrap();
            let point =
                |v: &ogeom_topo::Shape| model.node(v).unwrap().data().as_vertex().unwrap().point;
            Point::midpoint(point(&a), point(&b)).distance(at) < 1e-9
        })
        .expect("an edge there")
}

/// A block's three edges at one corner, rounded in one call at small
/// radii: the corner's ball patch and the bands meet along rims that leave
/// the block's sharp edges tangentially, and each closes to the closed
/// form: every band's section (1 - pi/4) r^2 over its length short of the
/// ball, and the ball's octant r^3 (1 - pi/6).
#[test]
fn a_corner_rounds_at_any_radius() {
    let (x, y, z) = (2.0, 2.0, 1.0);
    for radius in [0.01, 0.05, 0.1, 0.2] {
        let mut model = ogeom_topo::Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (x, y, z), T)
            .unwrap()
            .shape;
        let edges: Vec<ogeom_topo::Shape> = [
            Point::new(x, y / 2.0, z),
            Point::new(x, y, z / 2.0),
            Point::new(x / 2.0, y, z),
        ]
        .iter()
        .map(|&p| edge_through(&model, &block, p))
        .collect();
        let rounded = ogeom_fillet::fillet_edges(&mut model, &block, &edges, radius, T)
            .unwrap()
            .shape;
        assert!(ogeom_algo::check(&model, &rounded, T).unwrap().is_valid());
        let removed =
            (1.0 - core::f64::consts::FRAC_PI_4) * radius * radius * (x + y + z - 3.0 * radius)
                + radius.powi(3) * (1.0 - core::f64::consts::PI / 6.0);
        let measured = volume(&model, &rounded, 1e-5);
        let want = x * y * z - removed;
        assert!(
            (measured - want).abs() < 1e-6,
            "radius {radius}: {measured} against {want}"
        );
    }
}

/// Every edge of an L-bracket, in either order. The step's inner corner,
/// where its re-entrant edge meets two convex ones, is no convex corner,
/// and whichever blend comes first there takes it.
#[test]
fn every_edge_of_an_l_bracket_rounds_in_either_order() {
    for reversed in [false, true] {
        let mut model = ogeom_topo::Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T)
            .unwrap()
            .shape;
        let step = Frame::new(Point::new(1.0, -0.5, 1.0), Direction::Z, Direction::X, T).unwrap();
        let notch = ogeom_algo::make_box(&mut model, step, (2.0, 3.0, 2.0), T)
            .unwrap()
            .shape;
        let bracket = ogeom_bool::cut(&mut model, &block, &notch, T)
            .unwrap()
            .shape;
        let mut edges = explore_unique(&model, &bracket, ShapeType::Edge).unwrap();
        assert_eq!(edges.len(), 18);
        if reversed {
            edges.reverse();
        }
        let rounded = ogeom_fillet::fillet_edges(&mut model, &bracket, &edges, 0.1, T)
            .unwrap_or_else(|e| panic!("reversed {reversed}: {e}"))
            .shape;
        assert!(ogeom_algo::check(&model, &rounded, T).unwrap().is_valid());
        let v = volume(&model, &rounded, 1e-4);
        assert!((v - 6.0).abs() < 0.1, "reversed {reversed}: {v}");
    }
}

/// Every edge of a slab rounded at half its thickness: the blends on each
/// side meet along its middle and the side is gone, and each short edge
/// lies wholly inside the balls of the two corners at its ends, which
/// round it. What is left is the slab's mid-plane rectangle grown by the
/// radius, whose volume Steiner's formula gives: area times thickness,
/// a half-cylinder round the perimeter, and a ball.
#[test]
fn a_slab_rounds_fully_at_half_its_thickness() {
    let pi = core::f64::consts::PI;
    for (x, y) in [(2.0, 2.0), (3.0, 2.0)] {
        let mut model = ogeom_topo::Model::new();
        let slab = ogeom_algo::make_box(&mut model, Frame::WORLD, (x, y, 1.0), T)
            .unwrap()
            .shape;
        let edges = explore_unique(&model, &slab, ShapeType::Edge).unwrap();
        let rounded = ogeom_fillet::fillet_edges(&mut model, &slab, &edges, 0.5, T)
            .unwrap_or_else(|e| panic!("{x} x {y}: {e}"));
        let diagnosis = ogeom_algo::check(&model, &rounded.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{x} x {y}: {:?}", diagnosis.problems);
        let (a, b) = (x - 1.0, y - 1.0);
        let want = a * b + (a + b) * pi * 0.25 + pi / 6.0;
        let got = volume(&model, &rounded.shape, 1e-3);
        assert!((got - want).abs() < 1e-6, "{x} x {y}: {got} against {want}");
        for edge in &edges {
            assert!(rounded.history.is_deleted(edge));
        }
    }
}

/// The top rims of a plate's nine bores, rounded as one chain: no rim
/// meets another, so the blends are applied together. Each removes the
/// corner between the plate and its bore turned about the bore's axis:
/// area `r^2 (1 - pi/4)`, its centroid `r (10 - 3 pi) / (3 (4 - pi))` out
/// from the corner, by Pappus. Every rim is gone and its blend credited to
/// it.
#[test]
fn rims_that_meet_nothing_round_together() {
    let pi = core::f64::consts::PI;
    let mut model = ogeom_topo::Model::new();
    let plate = ogeom_algo::make_box(&mut model, Frame::WORLD, (15.0, 15.0, 4.0), T)
        .unwrap()
        .shape;
    let mut drums = Vec::new();
    for i in 0..3 {
        for j in 0..3 {
            let at = Frame::new(
                Point::new(2.5 + 5.0 * f64::from(i), 2.5 + 5.0 * f64::from(j), -1.0),
                Direction::Z,
                Direction::X,
                T,
            )
            .unwrap();
            drums.push(
                ogeom_algo::make_cylinder(&mut model, at, 1.0, 6.0, T)
                    .unwrap()
                    .shape,
            );
        }
    }
    let grid = model.add_compound(&drums).unwrap();
    let bored = ogeom_bool::cut(&mut model, &plate, &grid, T).unwrap().shape;
    let rims: Vec<_> = explore_unique(&model, &bored, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            use ogeom_geom::Curve3d as _;
            let (curve, range) = edge_curve_of(&model, e);
            matches!(curve, ogeom_geom::Curve::Circle(_))
                && curve
                    .point_at(f64::midpoint(range.0, range.1), T)
                    .unwrap()
                    .z
                    > 3.99
        })
        .collect();
    let r = 0.3;
    let rounded = ogeom_fillet::fillet_edges(&mut model, &bored, &rims, r, T).unwrap();
    assert!(
        ogeom_algo::check(&model, &rounded.shape, T)
            .unwrap()
            .is_valid()
    );
    let area = r * r * (1.0 - pi / 4.0);
    let centroid = 1.0 + r * (10.0 - 3.0 * pi) / (3.0 * (4.0 - pi));
    let removed = volume(&model, &bored, 1e-3) - volume(&model, &rounded.shape, 1e-3);
    let want = 9.0 * area * 2.0 * pi * centroid;
    assert!(
        (removed - want).abs() < want * 1e-4,
        "{removed} against {want}"
    );
    let mut credited = 0;
    for rim in &rims {
        assert!(rounded.history.is_deleted(rim));
        credited += rounded.history.generated(rim).len();
    }
    assert!(credited >= 9, "{credited}");
}

/// An edge's curve and range, read off the model.
fn edge_curve_of(
    model: &ogeom_topo::Model,
    edge: &ogeom_topo::Shape,
) -> (ogeom_geom::Curve, (f64, f64)) {
    let data = model.node(edge).unwrap().data().as_edge().unwrap();
    let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        panic!("an edge with no curve");
    };
    (model.geometry().curve(*curve).unwrap().clone(), *range)
}
