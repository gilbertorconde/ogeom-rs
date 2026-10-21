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

/// Every edge of a box standing at a slant, rounded: each corner's ball
/// meets its bands along meridians through the poles of the ball's own
/// chart, which the slant leaves a few nanoradians off the poles found
/// by projection. The volume is Steiner's: the box shrunk by the radius,
/// a slab on each face, a quarter cylinder along each edge and the ball
/// at the corners.
#[test]
fn every_edge_of_a_slanted_box_rounds() {
    let pi = core::f64::consts::PI;
    let z = Direction::new(Vector::new(0.3, -0.5, 0.8), T).unwrap();
    let x = Direction::new(z.vector().cross(Vector::new(0.0, 0.0, 1.0)), T).unwrap();
    let frame = Frame::new(Point::new(0.7, -1.1, 2.3), z, x, T).unwrap();
    for (dims, radius) in [
        ((2.0, 2.0, 2.0), 0.2),
        ((3.0, 2.0, 1.5), 0.15),
        ((4.0, 2.5, 3.3), 0.5),
    ] {
        let mut model = ogeom_topo::Model::new();
        let block = ogeom_algo::make_box(&mut model, frame, dims, T)
            .unwrap()
            .shape;
        let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
        let rounded = ogeom_fillet::fillet_edges(&mut model, &block, &edges, radius, T)
            .unwrap_or_else(|e| panic!("{dims:?} at {radius}: {e}"));
        let diagnosis = ogeom_algo::check(&model, &rounded.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{dims:?}: {:?}", diagnosis.problems);
        let (a, b, c) = (
            dims.0 - 2.0 * radius,
            dims.1 - 2.0 * radius,
            dims.2 - 2.0 * radius,
        );
        let want = a * b * c
            + 2.0 * radius * (a * b + b * c + a * c)
            + pi * radius * radius * (a + b + c)
            + 4.0 / 3.0 * pi * radius.powi(3);
        let got = volume(&model, &rounded.shape, 1e-4);
        assert!(
            (got - want).abs() < 1e-5 * want,
            "{dims:?}: {got} against {want}"
        );
    }
}

/// A slanted box whose corner tool cuts its block's planes along lines
/// running the full unbounded domain: a step toward a line's end lands an
/// ulp of a billion past it, which is the end, not beyond it.
#[test]
fn every_edge_of_a_box_whose_sections_run_to_a_line_s_end_rounds() {
    let pi = core::f64::consts::PI;
    let z = Direction::new(
        Vector::new(0.5622960444626135, 0.5069790817054428, 0.6532957745882064),
        T,
    )
    .unwrap();
    let x = Direction::new(Vector::new(0.0, 0.7900197472052799, -0.6130813967375831), T).unwrap();
    let at = Point::new(0.9018390493041868, 1.3163269939520887, 1.694201457746387);
    let frame = Frame::new(at, z, x, T).unwrap();
    let dims = (1.8161039209587628, 2.908987621268594, 3.5473365999764437);
    let radius = 0.19069391608246453;
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, frame, dims, T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rounded = ogeom_fillet::fillet_edges(&mut model, &block, &edges, radius, T)
        .unwrap_or_else(|e| panic!("{e}"));
    let diagnosis = ogeom_algo::check(&model, &rounded.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let (a, b, c) = (
        dims.0 - 2.0 * radius,
        dims.1 - 2.0 * radius,
        dims.2 - 2.0 * radius,
    );
    let want = a * b * c
        + 2.0 * radius * (a * b + b * c + a * c)
        + pi * radius * radius * (a + b + c)
        + 4.0 / 3.0 * pi * radius.powi(3);
    let got = volume(&model, &rounded.shape, 1e-4);
    assert!((got - want).abs() < 1e-5 * want, "{got} against {want}");
}

/// A box a hair off the world axes, every edge rounded. The ball's rims
/// are circles the boolean starts where a world axis falls in their plane,
/// here three ten-thousandths of a radian from the touch points; split
/// there, the sliver of arc between was taken for dust and the arc beside
/// it left short of its neighbours. The volume is Steiner's.
#[test]
fn every_edge_of_a_box_a_hair_off_the_axes_rounds() {
    let pi = core::f64::consts::PI;
    let z = Direction::new(
        Vector::new(
            0.0003048220105320179,
            0.9836284277769499,
            0.18020827714782525,
        ),
        T,
    )
    .unwrap();
    let x = Direction::new(
        Vector::new(
            -0.9999999535417698,
            0.0002998316089010636,
            5.4931451906729234e-5,
        ),
        T,
    )
    .unwrap();
    let at = Point::new(
        -0.09450952519042044,
        -1.411639443290692,
        -1.0075285054870318,
    );
    let frame = Frame::new(at, z, x, T).unwrap();
    let dims = (2.03115597451817, 3.7040754186675153, 1.5324456764276202);
    let radius = 0.4739374561910516;
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, frame, dims, T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rounded = ogeom_fillet::fillet_edges(&mut model, &block, &edges, radius, T)
        .unwrap_or_else(|e| panic!("{e}"));
    let diagnosis = ogeom_algo::check(&model, &rounded.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let (a, b, c) = (
        dims.0 - 2.0 * radius,
        dims.1 - 2.0 * radius,
        dims.2 - 2.0 * radius,
    );
    let want = a * b * c
        + 2.0 * radius * (a * b + b * c + a * c)
        + pi * radius * radius * (a + b + c)
        + 4.0 / 3.0 * pi * radius.powi(3);
    let got = volume(&model, &rounded.shape, 1e-4);
    assert!((got - want).abs() < 1e-5 * want, "{got} against {want}");
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

/// A keyhole: a drum of radius 2, its seam at 310 degrees, fused with a
/// neck 2 wide, 1 high. The drum's wall is one face, and its rim on the
/// top splits at the seam into a long arc and a short one, 310 to 330
/// degrees. The rim stops where the neck leaves it, so the short arc's
/// round runs its own 20 degrees and no further, though the wall and the
/// long arc stand opposite its middle. Rounded alone, and with the long
/// arc, each round removes the corner between the top and the wall
/// turned about the axis over its arc: area `r^2 (1 - pi/4)`, its
/// centroid `r (10 - 3 pi) / (3 (4 - pi))` in from the rim, by Pappus.
#[test]
fn a_short_arc_of_a_rim_that_stops_rounds_over_its_own_arc() {
    let pi = core::f64::consts::PI;
    let radius = 2.0;
    let mut model = ogeom_topo::Model::new();
    let at = |degrees: f64| {
        let a = degrees.to_radians();
        Point::new(radius * a.cos(), radius * a.sin(), 0.0)
    };
    let seam = Direction::new(at(310.0) - Point::ORIGIN, T).unwrap();
    let drum = ogeom_algo::make_cylinder(
        &mut model,
        Frame::new(Point::ORIGIN, Direction::Z, seam, T).unwrap(),
        radius,
        1.0,
        T,
    )
    .unwrap()
    .shape;
    let neck = ogeom_algo::make_box(
        &mut model,
        Frame::new(Point::new(0.0, -1.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
        (5.0, 2.0, 1.0),
        T,
    )
    .unwrap()
    .shape;
    let keyhole = ogeom_bool::fuse(&mut model, &drum, &neck, T).unwrap().shape;
    // The top edge whose bounds hold `inside` and lie within `within` of
    // it: the long arc alone reaches round to x = -2, and the short one
    // alone stays near 320 degrees.
    let top = |inside: Point, within: f64| {
        explore_unique(&model, &keyhole, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let bounds = ogeom_algo::shape_bounds(&model, e, T).unwrap();
                let (lo, hi) = (bounds.low().unwrap(), bounds.high().unwrap());
                lo.z > 0.999
                    && bounds.distance_to(inside) < 1e-6
                    && lo.distance(inside).max(hi.distance(inside)) < within
            })
            .expect("the rim has that arc")
    };
    let middle = at(320.0);
    let (long, short) = (
        top(Point::new(-radius, 0.0, 1.0), 5.0),
        top(Point::new(middle.x, middle.y, 1.0), 0.5),
    );
    let r = 0.3;
    let area = r * r * (1.0 - pi / 4.0);
    let centroid = radius - r * (10.0 - 3.0 * pi) / (3.0 * (4.0 - pi));
    let before = volume(&model, &keyhole, 1e-4);
    for (named, sweep) in [(vec![short.clone()], 20.0), (vec![long, short], 300.0)] {
        let mut model = model.clone();
        let rounded = ogeom_fillet::fillet_edges(&mut model, &keyhole, &named, r, T).unwrap();
        assert!(
            ogeom_algo::check(&model, &rounded.shape, T)
                .unwrap()
                .is_valid()
        );
        let removed = before - volume(&model, &rounded.shape, 1e-4);
        let want = area * centroid * f64::to_radians(sweep);
        assert!(
            (removed - want).abs() < want * 1e-3,
            "{sweep} degrees: {removed} against {want}"
        );
    }
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

/// What two rounds of `radius` meeting across a square inside corner take
/// past their caps: each runs on to the plane bisecting the corner, and the
/// section a distance w in from its wall reaches w along, so each adds the
/// corner outside its quarter circle weighted by w, r^3 (5/6 - pi/4).
fn square_corner_cell(radius: f64) -> f64 {
    2.0 * (5.0 / 6.0 - core::f64::consts::FRAC_PI_4) * radius.powi(3)
}

/// An L-shaped block whose top has an inside corner at `(10, 10, 5)`, and
/// its two top edges leaving that corner, the first along `y = 10`. The
/// first edge's curve runs `overrun` past the corner, and the corner's
/// vertex tolerance is widened to `reach` to cover it. Where
/// `vertex_follows`, the vertex stands where the stretched curve ends, off
/// the wall `x = 10` the first edge meets there and off the second edge's
/// line, rather than at the corner.
fn l_block_with_an_edge_past_its_inside_corner(
    model: &mut ogeom_topo::Model,
    overrun: f64,
    reach: f64,
    vertex_follows: bool,
) -> (ogeom_topo::Shape, [ogeom_topo::Shape; 2]) {
    let block = ogeom_algo::make_box(model, Frame::WORLD, (20.0, 20.0, 5.0), T)
        .unwrap()
        .shape;
    let notch_at = Frame::new(Point::new(10.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let notch = ogeom_algo::make_box(model, notch_at, (20.0, 20.0, 7.0), T)
        .unwrap()
        .shape;
    let l_block = ogeom_bool::cut(model, &block, &notch, T).unwrap().shape;
    let edges = [
        edge_through(model, &l_block, Point::new(15.0, 10.0, 5.0)),
        edge_through(model, &l_block, Point::new(10.0, 15.0, 5.0)),
    ];
    let corner = Point::new(10.0, 10.0, 5.0);
    // The parameter of the first edge's curve at the corner, and how far
    // past it the curve's range is stretched: the edge's representations
    // share one parameter, so each range moves at that end alike.
    let (curve, (lo, hi)) = edge_curve_of(model, &edges[0]);
    let at = |t: f64| ogeom_geom::Curve3d::point_at(&curve, t, T).unwrap();
    let speed = at(lo).distance(at(hi)) / (hi - lo);
    let (end, past) = if at(lo).distance(corner) < at(hi).distance(corner) {
        (lo, lo - overrun / speed)
    } else {
        (hi, hi + overrun / speed)
    };
    let node = model.node_mut(&edges[0]).unwrap();
    let ogeom_topo::NodeData::Edge(data) = node.data_mut() else {
        panic!("an edge node holds edge data");
    };
    for repr in &mut data.representations {
        if let ogeom_topo::EdgeRepr::Curve3d { range, .. }
        | ogeom_topo::EdgeRepr::PCurve { range, .. }
        | ogeom_topo::EdgeRepr::Seam { range, .. } = repr
        {
            for bound in [&mut range.0, &mut range.1] {
                if (*bound - end).abs() < 1e-12 {
                    *bound = past;
                }
            }
        }
    }
    let vertex = explore_unique(model, &l_block, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .find(|v| {
            model
                .node(v)
                .unwrap()
                .data()
                .as_vertex()
                .unwrap()
                .point
                .distance(corner)
                < 1e-9
        })
        .expect("the inside corner's vertex");
    model
        .widen(&vertex, ogeom_core::Tolerance::new(reach).unwrap())
        .unwrap();
    if vertex_follows {
        let ogeom_topo::NodeData::Vertex(data) = model.node_mut(&vertex).unwrap().data_mut() else {
            panic!("a vertex node holds vertex data");
        };
        data.point = at(past);
    }
    (l_block, edges)
}

/// Two top edges of an L-shaped block rounded, and then chamfered, where
/// they meet at the top's inside corner, one edge's curve running a little
/// past the corner's vertex, within the vertex's tolerance. Each blend
/// removes its edge's section over the edge's length from the outside wall
/// to the inside corner, in either order; the two rounds also meet in the
/// plane bisecting the corner, each running on past its cap to it, which
/// takes the cell past both caps that both hold, two times r^3 (5/6 -
/// pi/4). The bevels stand capped at the vertex, flush with each other. A
/// cap left where the curve ends would cut into the other blend past the
/// corner.
#[test]
fn blends_at_an_inside_corner_their_edge_runs_past_meet_at_the_vertex() {
    let size = 1.0;
    // Far above the confusion, so a cap left where the curve ends would
    // cut a measurable notch, and far below the blends' size.
    let overrun = 1e-3;
    let round = (1.0 - core::f64::consts::FRAC_PI_4) * size * size;
    let bevel = size * size / 2.0;
    for (chamfer, section) in [(false, round), (true, bevel)] {
        for reversed in [false, true] {
            let mut model = ogeom_topo::Model::new();
            let (block, mut edges) = l_block_with_an_edge_past_its_inside_corner(
                &mut model,
                overrun,
                overrun * 2.0,
                false,
            );
            if reversed {
                edges.reverse();
            }
            let blended = if chamfer {
                ogeom_fillet::chamfer_edges(&mut model, &block, &edges, size, T)
            } else {
                ogeom_fillet::fillet_edges(&mut model, &block, &edges, size, T)
            }
            .unwrap_or_else(|e| panic!("chamfer {chamfer}, reversed {reversed}: {e}"))
            .shape;
            assert!(ogeom_algo::check(&model, &blended, T).unwrap().is_valid());
            let cell = if chamfer {
                0.0
            } else {
                square_corner_cell(size)
            };
            let want = 1500.0 - section * 20.0 - cell;
            let measured = volume(&model, &blended, 1e-5);
            // The mesh at a hundred-thousandth chord reads a box corner's
            // rounds to a micro-cubic-millimetre
            // (a_corner_rounds_at_any_radius), and a bevel's planes exactly;
            // the notch an overrun cap cuts is its section times 1e-3.
            assert!(
                (measured - want).abs() < 2e-5,
                "chamfer {chamfer}, reversed {reversed}: {measured} against {want}"
            );
        }
    }
}

/// The same two blends where the corner's vertex stands at the end of the
/// first edge's curve, which runs past the wall `x = 10` that ends it, and
/// off the second edge's line by as much. Measured from the vertex, neither
/// line runs past it; measured against the second edge, which runs along
/// that wall, the first does. The blends meet where the edges' lines do:
/// the bevels capped there, flush with each other, and the rounds in the
/// plane bisecting the corner through it, as above, in either order. A cap
/// at the vertex would cut into the other blend past the corner.
#[test]
fn blends_at_an_inside_corner_whose_vertex_stands_off_the_wall_meet_at_the_wall() {
    let size = 1.0;
    // As in the test above: a notch of the section times the overrun is
    // far above what the mesh reads, and the overrun far below the blends.
    let overrun = 1e-3;
    let round = (1.0 - core::f64::consts::FRAC_PI_4) * size * size;
    let bevel = size * size / 2.0;
    for (chamfer, section) in [(false, round), (true, bevel)] {
        for reversed in [false, true] {
            let mut model = ogeom_topo::Model::new();
            let (block, mut edges) = l_block_with_an_edge_past_its_inside_corner(
                &mut model,
                overrun,
                overrun * 2.0,
                true,
            );
            if reversed {
                edges.reverse();
            }
            let blended = if chamfer {
                ogeom_fillet::chamfer_edges(&mut model, &block, &edges, size, T)
            } else {
                ogeom_fillet::fillet_edges(&mut model, &block, &edges, size, T)
            }
            .unwrap_or_else(|e| panic!("chamfer {chamfer}, reversed {reversed}: {e}"))
            .shape;
            assert!(ogeom_algo::check(&model, &blended, T).unwrap().is_valid());
            let cell = if chamfer {
                0.0
            } else {
                square_corner_cell(size)
            };
            let want = 1500.0 - section * 20.0 - cell;
            let measured = volume(&model, &blended, 1e-5);
            // The same reading as the test above: the notch a cap past the
            // wall cuts is the section times 1e-3, ten times this bound.
            assert!(
                (measured - want).abs() < 2e-5,
                "chamfer {chamfer}, reversed {reversed}: {measured} against {want}"
            );
        }
    }
}

/// A prism two high whose top has a re-entrant corner at `(4, 3)`: its
/// outline runs along `y = 3` from `x = 8` to the corner and turns `turn`
/// back towards the material there for three along, each of the two edges
/// ending square on an end wall. Returns the solid and the two top edges
/// leaving the corner.
fn re_entrant_corner(
    model: &mut ogeom_topo::Model,
    turn: f64,
) -> (ogeom_topo::Shape, [ogeom_topo::Shape; 2]) {
    let (s, c) = turn.sin_cos();
    let along = Vector::new(-c, s, 0.0);
    let corner = Point::new(4.0, 3.0, 0.0);
    let end = corner + along * 3.0;
    let foot = end + Vector::new(-s, -c, 0.0);
    let outline = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(8.0, 0.0, 0.0),
        Point::new(8.0, 3.0, 0.0),
        corner,
        end,
        foot,
        Point::new(0.0, foot.y, 0.0),
    ];
    let wire = ogeom_algo::make_polygon(model, &outline, true, T)
        .unwrap()
        .shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let prism = ogeom_algo::make_prism(model, &face, Vector::new(0.0, 0.0, 2.0), T)
        .unwrap()
        .shape;
    let placed_edge_through = |at: Point| {
        explore_unique(model, &prism, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let (a, b) = ogeom_algo::edge_vertices(model, e).unwrap().unwrap();
                let point = |v: &ogeom_topo::Shape| {
                    let p = model.node(v).unwrap().data().as_vertex().unwrap().point;
                    v.transform(model.datums()).unwrap().apply(p)
                };
                Point::midpoint(point(&a), point(&b)).distance(at) < 1e-9
            })
            .expect("a top edge there")
    };
    let edges = [
        placed_edge_through(Point::new(6.0, 3.0, 2.0)),
        placed_edge_through(corner + along * 1.5 + Vector::new(0.0, 0.0, 2.0)),
    ];
    (prism, edges)
}

/// Two top edges rounded where the outline turns back into the material
/// between them, at a square corner and at a shallower one: the rounds
/// meet in the plane bisecting the corner, each running on past its cap
/// to it, so the section a distance w in from its wall reaches w tan(t/2)
/// along, t the turn. Together they take each edge's section over its
/// length, the corner outside the quarter circle, and twice that corner's
/// moment about its wall times tan(t/2), in either order.
#[test]
fn rounds_across_a_re_entrant_turn_meet_in_the_plane_bisecting_it() {
    let pi = core::f64::consts::PI;
    let r = 0.5;
    let section = (1.0 - pi / 4.0) * r * r;
    let moment = (5.0 / 6.0 - pi / 4.0) * r.powi(3);
    for degrees in [90.0_f64, 45.0, 150.0] {
        let turn = degrees.to_radians();
        for reversed in [false, true] {
            let mut model = ogeom_topo::Model::new();
            let (prism, mut edges) = re_entrant_corner(&mut model, turn);
            let before = volume(&model, &prism, 1e-5);
            if reversed {
                edges.reverse();
            }
            let rounded = ogeom_fillet::fillet_edges(&mut model, &prism, &edges, r, T)
                .unwrap_or_else(|e| panic!("{degrees} degrees, reversed {reversed}: {e}"))
                .shape;
            assert!(ogeom_algo::check(&model, &rounded, T).unwrap().is_valid());
            let removed = before - volume(&model, &rounded, 1e-5);
            let want = section.mul_add(7.0, 2.0 * moment * (turn / 2.0).tan());
            assert!(
                (removed - want).abs() < 2e-5,
                "{degrees} degrees, reversed {reversed}: {removed} against {want}"
            );
        }
    }
}

/// The volume `fillet_edges` takes off `solid` rounding `edges` at
/// `radius`, in the order given and reversed, the result checked.
fn removed_both_ways(
    model: &ogeom_topo::Model,
    solid: &ogeom_topo::Shape,
    edges: &[ogeom_topo::Shape],
    radius: f64,
) -> [f64; 2] {
    let before = volume(model, solid, 1e-4);
    let mut reversed = edges.to_vec();
    reversed.reverse();
    [edges.to_vec(), reversed].map(|named| {
        let mut model = model.clone();
        let rounded = ogeom_fillet::fillet_edges(&mut model, solid, &named, radius, T)
            .unwrap()
            .shape;
        assert!(ogeom_algo::check(&model, &rounded, T).unwrap().is_valid());
        before - volume(&model, &rounded, 1e-4)
    })
}

/// The edge of `shape` whose ends' midpoint stands at `at`, its vertices
/// read where they are placed.
fn placed_edge_through(
    model: &ogeom_topo::Model,
    shape: &ogeom_topo::Shape,
    at: Point,
) -> ogeom_topo::Shape {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (a, b) = ogeom_algo::edge_vertices(model, e).unwrap().unwrap();
            let point = |v: &ogeom_topo::Shape| {
                let p = model.node(v).unwrap().data().as_vertex().unwrap().point;
                v.transform(model.datums()).unwrap().apply(p)
            };
            Point::midpoint(point(&a), point(&b)).distance(at) < 1e-6
        })
        .expect("an edge there")
}

/// A face on the plane through `on` with unit `normal` bounded by the
/// polygon `outline`, swept along `along`.
fn polygon_prism(
    model: &mut ogeom_topo::Model,
    outline: &[Point],
    normal: Direction,
    along: Vector,
) -> ogeom_topo::Shape {
    let wire = ogeom_algo::make_polygon(model, outline, true, T)
        .unwrap()
        .shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(outline[0], normal),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    ogeom_algo::make_prism(model, &face, along, T)
        .unwrap()
        .shape
}

/// A round's section across a straight crease between a host `z <= 2`
/// and a wall through the crease at `(edge, 2)` in a plane `(u, z)` square
/// to it, the wall's outward normal `(nu, nz)`: the material lies under the
/// host and behind the wall. Measured a column at a time: the heights `z`
/// a vertical line at `u` holds of the round's section, and of what stands
/// past the wall's plane under the host.
struct Section {
    centre: (f64, f64),
    radius: f64,
    /// Half-planes `a.(u, z) <= b` bounding the corner the ball rounds.
    sides: Vec<((f64, f64), f64)>,
    /// The wall's half-plane, the material on its `<=` side.
    wall: ((f64, f64), f64),
    area: f64,
}

impl Section {
    fn new(edge: f64, (nu, nz): (f64, f64), radius: f64) -> Self {
        let m = nu.hypot(nz);
        let n = (nu / m, nz / m);
        let crease = (edge, 2.0);
        // The ball's centre a radius in from both the host and the wall.
        let k = radius / (1.0 + n.1);
        let centre = (crease.0 - n.0 * k, crease.1 - (1.0 + n.1) * k);
        let wall = (n, n.0.mul_add(crease.0, n.1 * crease.1));
        let mut sides = vec![((0.0, 1.0), 2.0), wall];
        // Each contact's radius bounds the corner on the crease's side.
        for (cu, cz) in [(0.0, 1.0), n] {
            let (mut au, mut az) = (cz, -cu);
            if au.mul_add(crease.0 - centre.0, az * (crease.1 - centre.1)) < 0.0 {
                (au, az) = (-au, -az);
            }
            sides.push(((-au, -az), -au.mul_add(centre.0, az * centre.1)));
        }
        let dihedral = (-n.1).acos();
        let area = radius * radius / (dihedral / 2.0).tan()
            - radius * radius * (core::f64::consts::PI - dihedral) / 2.0;
        Self {
            centre,
            radius,
            sides,
            wall,
            area,
        }
    }

    /// The heights of the column at `u` in the section: under the host and
    /// within the corner, less the ball.
    fn column(&self, u: f64) -> [(f64, f64); 2] {
        let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
        for &((au, az), b) in &self.sides {
            if az.abs() < 1e-15 {
                if au * u > b {
                    return [(0.0, 0.0); 2];
                }
            } else if az > 0.0 {
                hi = hi.min((b - au * u) / az);
            } else {
                lo = lo.max((b - au * u) / az);
            }
        }
        if hi <= lo {
            return [(0.0, 0.0); 2];
        }
        let reach = self
            .radius
            .mul_add(self.radius, -(u - self.centre.0).powi(2));
        if reach <= 0.0 {
            return [(lo, hi), (0.0, 0.0)];
        }
        let h = reach.sqrt();
        [
            (lo, hi.min(self.centre.1 - h)),
            (lo.max(self.centre.1 + h), hi),
        ]
    }

    /// The heights of the column at `u` past the wall's plane, under the
    /// host, down to `floor`.
    fn past(&self, u: f64, floor: f64) -> (f64, f64) {
        let ((au, az), b) = self.wall;
        if az.abs() < 1e-15 {
            return if au * u > b { (floor, 2.0) } else { (0.0, 0.0) };
        }
        let at = (b - au * u) / az;
        if az > 0.0 {
            (at.max(floor), 2.0)
        } else {
            (floor, at.min(2.0))
        }
    }
}

/// How much of the heights `a` and `b`, each a few disjoint runs, overlap.
fn overlap(a: &[(f64, f64)], b: &[(f64, f64)]) -> f64 {
    let mut sum = 0.0;
    for &(p, q) in a {
        for &(r, t) in b {
            sum += (q.min(t) - p.max(r)).max(0.0);
        }
    }
    sum
}

/// The total length of the heights `a`, a few disjoint runs.
fn span(a: &[(f64, f64)]) -> f64 {
    a.iter().map(|&(p, q)| (q - p).max(0.0)).sum()
}

/// The midpoint sum of `f` over the rectangle `[x0, x1] x [y0, y1]`, `n`
/// cells to a unit each way.
fn integral(
    (x0, x1): (f64, f64),
    (y0, y1): (f64, f64),
    n: f64,
    f: impl Fn(f64, f64) -> f64,
) -> f64 {
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "a cell count"
    )]
    let (nx, ny) = (
        ((x1 - x0) * n).round() as usize,
        ((y1 - y0) * n).round() as usize,
    );
    #[allow(clippy::cast_precision_loss, reason = "a cell count")]
    let (dx, dy) = ((x1 - x0) / nx as f64, (y1 - y0) / ny as f64);
    let mut sum = 0.0;
    for i in 0..nx {
        #[allow(clippy::cast_precision_loss, reason = "a cell index")]
        let x = (i as f64 + 0.5).mul_add(dx, x0);
        for j in 0..ny {
            #[allow(clippy::cast_precision_loss, reason = "a cell index")]
            let y = (j as f64 + 0.5).mul_add(dy, y0);
            sum += f(x, y);
        }
    }
    sum * dx * dy
}

/// What the corner of [`drafted_inside_corner`] takes beyond its two
/// creases' sections along their lengths, its drafted wall's section and
/// its square one's given: past both caps what lies above both bands, and
/// before the square wall's cap what of its section lies above the drafted
/// band less that section, a column at a time.
fn drafted_corner_gain(drafted: &Section, square: &Section, draft: f64) -> f64 {
    integral((2.0, 3.0), (2.0, 3.0 + draft.abs()), 2e3, |x, y| {
        let [a0, a1] = drafted.column(y);
        let [b0, b1] = square.column(x);
        let above = [a0, a1, drafted.past(y, -1.0)];
        if y < 3.0 {
            overlap(&above, &[b0, b1, square.past(x, -1.0)])
        } else {
            overlap(&above, &[b0, b1]) - span(&[b0, b1])
        }
    })
}

/// A block six square and two high with a notch where `x >= 3` and `y` is
/// past the wall `y = 3 - d + d z / 2`, drafted by `d` over its height,
/// all turned `spin` about the `z` axis: returns the block and the two top
/// creases leaving the notch's inside corner `(3, 3, 2)` (turned), the
/// drafted wall's first.
fn drafted_inside_corner(
    model: &mut ogeom_topo::Model,
    draft: f64,
    spin: f64,
) -> (ogeom_topo::Shape, [ogeom_topo::Shape; 2]) {
    let (s, c) = spin.sin_cos();
    let turn = |p: Point| Point::new(p.x * c - p.y * s, p.x * s + p.y * c, p.z);
    let turn_v = |v: Vector| Vector::new(v.x * c - v.y * s, v.x * s + v.y * c, v.z);
    let frame = Frame::new(
        Point::ORIGIN,
        Direction::Z,
        Direction::new(turn_v(Vector::X), T).unwrap(),
        T,
    )
    .unwrap();
    let block = ogeom_algo::make_box(model, frame, (6.0, 6.0, 2.0), T)
        .unwrap()
        .shape;
    let wall = |z: f64| draft.mul_add(z / 2.0, 3.0 - draft);
    let notch = polygon_prism(
        model,
        &[
            turn(Point::new(3.0, wall(-1.0), -1.0)),
            turn(Point::new(3.0, 7.0, -1.0)),
            turn(Point::new(3.0, 7.0, 3.0)),
            turn(Point::new(3.0, wall(3.0), 3.0)),
        ],
        Direction::new(turn_v(Vector::X), T).unwrap(),
        turn_v(Vector::new(4.0, 0.0, 0.0)),
    );
    let cut = ogeom_bool::cut(model, &block, &notch, T).unwrap().shape;
    let block = explore_unique(model, &cut, ShapeType::Solid).unwrap()[0].clone();
    let edges = [
        placed_edge_through(model, &block, turn(Point::new(4.5, 3.0, 2.0))),
        placed_edge_through(model, &block, turn(Point::new(3.0, 4.5, 2.0))),
    ];
    (block, edges)
}

/// Whether every face `edge` generated lies along it: on a cylinder or an
/// extrusion whose direction is the edge's own (its band), or on a plane
/// holding that direction (a host's plane carried past it).
fn generated_along(
    model: &ogeom_topo::Model,
    history: &ogeom_algo::History,
    edge: &ogeom_topo::Shape,
    along: Vector,
) -> bool {
    let faces = history.generated(edge);
    !faces.is_empty()
        && faces.iter().all(|face| {
            let data = model.node(face).unwrap().data().as_face().unwrap();
            let placed = face.transform(model.datums()).unwrap();
            match model.geometry().surface(data.surface).unwrap() {
                ogeom_geom::SurfaceGeometry::Cylinder(c) => {
                    placed
                        .apply_vector(c.cylinder().axis().direction.vector())
                        .cross(along)
                        .magnitude()
                        < 1e-9
                }
                ogeom_geom::SurfaceGeometry::Extrusion(e) => {
                    placed
                        .apply_vector(e.direction().vector())
                        .cross(along)
                        .magnitude()
                        < 1e-9
                }
                ogeom_geom::SurfaceGeometry::Plane(p) => {
                    placed
                        .apply_vector(p.plane().normal().vector())
                        .dot(along)
                        .abs()
                        < 1e-9
                }
                _ => false,
            }
        })
}

/// Two top edges of a block rounded where they meet at its inside corner,
/// one wall there drafted so the two bands' end sections differ, the
/// overhanging way and the moulded way. Each band carried on is the
/// straight band of its end section, above which stands what its ball
/// rolls off and what lies past its wall's plane; the corner takes what
/// lies above both. Measured against that region a column at a time: each
/// crease's section along it, and near the corner, past both caps what
/// stands above both bands, and before the square wall's cap what of its
/// section stands above the drafted band. Capped at the vertex instead,
/// the bands left the corner cell and, the moulded way, the square band
/// cut under the drafted one. Each band's faces stay its own edge's.
#[test]
fn rounds_at_a_drafted_inside_corner_meet_along_their_own_intersection() {
    let r = 0.5;
    let spun = 30.0_f64.to_radians();
    for (draft, spin) in [
        (0.005_f64, 0.0),
        (0.01, 0.0),
        (0.02, 0.0),
        (0.05, 0.0),
        (0.2, 0.0),
        (-0.2, 0.0),
        (-0.5, 0.0),
        (0.2, spun),
        (-0.2, spun),
    ] {
        let drafted = Section::new(3.0, (1.0, -draft / 2.0), r);
        let square = Section::new(3.0, (1.0, 0.0), r);
        let near = drafted_corner_gain(&drafted, &square, draft);
        let want = 3.0f64.mul_add(drafted.area + square.area, near);
        let mut model = ogeom_topo::Model::new();
        let (block, edges) = drafted_inside_corner(&mut model, draft, spin);
        let before = volume(&model, &block, 1e-4);
        for named in [edges.clone(), [edges[1].clone(), edges[0].clone()]] {
            let mut model = model.clone();
            let rounded = ogeom_fillet::fillet_edges(&mut model, &block, &named, r, T).unwrap();
            assert!(
                ogeom_algo::check(&model, &rounded.shape, T)
                    .unwrap()
                    .is_valid()
            );
            let removed = before - volume(&model, &rounded.shape, 1e-4);
            assert!(
                (removed - want).abs() < 3e-6,
                "draft {draft} spin {spin}: {removed} against {want}"
            );
            let (s, c) = spin.sin_cos();
            for (edge, along) in edges
                .iter()
                .zip([Vector::new(c, s, 0.0), Vector::new(-s, c, 0.0)])
            {
                assert!(
                    generated_along(&model, &rounded.history, edge, along),
                    "draft {draft} spin {spin}: an edge credited a face off its own band"
                );
            }
        }
    }
}

/// Two top edges of a block rounded where they meet at its outside corner,
/// the side wall drafted so its top stands `d` in. The two end sections
/// differ, and each band runs out past the vertex as it does alone, along
/// its own straight continuation until the solid has left its section; the
/// two cuts take their union, in either order. Measured against that union
/// a column at a time near the corner. Built after the side wall's band,
/// the front band stopped at its vertex where the solid ran on under it.
#[test]
fn rounds_at_a_drafted_outside_corner_take_both_bands_run_out() {
    let r = 0.5;
    for draft in [0.2_f64, 0.5, 0.8] {
        let front = Section::new(0.0, (1.0, 0.0), r);
        let side = Section::new(4.0 - draft, (1.0, draft / 2.0), r);
        let near = integral((2.6, 4.0), (0.0, 0.6), 2e3, |x, y| {
            // The solid's column: under the top and behind the side wall.
            let solid = (0.0, ((4.0 - x) * 2.0 / draft).min(2.0));
            let clip = |runs: [(f64, f64); 2]| runs.map(|(p, q)| (p.max(solid.0), q.min(solid.1)));
            let a = clip(front.column(-y));
            let b = clip(side.column(x));
            span(&a) + span(&b) - overlap(&a, &b)
        });
        let want = 2.6f64.mul_add(front.area, 3.4 * side.area) + near;
        let mut model = ogeom_topo::Model::new();
        let block = polygon_prism(
            &mut model,
            &[
                Point::new(0.0, 0.0, 0.0),
                Point::new(4.0, 0.0, 0.0),
                Point::new(4.0 - draft, 0.0, 2.0),
                Point::new(0.0, 0.0, 2.0),
            ],
            Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
            Vector::new(0.0, 4.0, 0.0),
        );
        let edges = [
            placed_edge_through(&model, &block, Point::new((4.0 - draft) / 2.0, 0.0, 2.0)),
            placed_edge_through(&model, &block, Point::new(4.0 - draft, 2.0, 2.0)),
        ];
        for got in removed_both_ways(&model, &block, &edges, r) {
            assert!(
                (got - want).abs() < 3e-6,
                "draft {draft}: {got} against {want}"
            );
        }
    }
}

/// A block eight square and two high with a pocket `[2, 6]` square and one
/// deep, its wall `y = 2` drafted by `d` over the pocket's depth: returns
/// the block and the pocket's four rim creases, the drafted wall's first,
/// then round the rim.
fn drafted_pocket(
    model: &mut ogeom_topo::Model,
    draft: f64,
) -> (ogeom_topo::Shape, Vec<ogeom_topo::Shape>) {
    let block = ogeom_algo::make_box(model, Frame::WORLD, (8.0, 8.0, 2.0), T)
        .unwrap()
        .shape;
    let wall = |z: f64| draft.mul_add(2.0 - z, 2.0);
    let pocket = polygon_prism(
        model,
        &[
            Point::new(2.0, wall(1.0), 1.0),
            Point::new(2.0, 6.0, 1.0),
            Point::new(2.0, 6.0, 3.0),
            Point::new(2.0, wall(3.0), 3.0),
        ],
        Direction::X,
        Vector::new(4.0, 0.0, 0.0),
    );
    let cut = ogeom_bool::cut(model, &block, &pocket, T).unwrap().shape;
    let block = explore_unique(model, &cut, ShapeType::Solid).unwrap()[0].clone();
    let edges = [(4.0, 2.0), (6.0, 4.0), (4.0, 6.0), (2.0, 4.0)]
        .map(|(x, y)| placed_edge_through(model, &block, Point::new(x, y, 2.0)))
        .to_vec();
    (block, edges)
}

/// A pocket's four rim creases rounded together, one wall drafted: the two
/// corners at the square walls meet flush in their bisecting planes, and
/// the two at the drafted wall are crossed in the same loop, the drafted
/// band cut in two halfway along and each half joined into its own corner.
/// Each crease takes its section along its length, each square corner the
/// cell of its two quarter rounds, 2 r^3 (5/6 - pi/4), and each drafted
/// corner what the drafted inside corner of the block above takes beyond
/// its sections, its draft over its height twice the pocket's per unit.
/// Capped at the drafted corners instead, the request left them standing.
/// With the pocket's drafted wall overhanging, the same holds.
#[test]
fn a_drafted_pocket_s_rim_meets_at_every_corner() {
    let r: f64 = 0.5;
    let cell = 2.0 * (5.0 / 6.0 - core::f64::consts::PI / 4.0) * r.powi(3);
    for draft in [0.1_f64, -0.1] {
        let mut model = ogeom_topo::Model::new();
        let (block, edges) = drafted_pocket(&mut model, draft);
        // The pocket's wall leans `draft` per unit of depth; the block's
        // drafted wall `-draft` per unit too, over its height of two.
        let block_draft = -2.0 * draft;
        let drafted = Section::new(3.0, (1.0, -block_draft / 2.0), r);
        let square = Section::new(3.0, (1.0, 0.0), r);
        let gain = drafted_corner_gain(&drafted, &square, block_draft);
        let want = 4.0f64.mul_add(
            3.0f64.mul_add(square.area, drafted.area),
            2.0f64.mul_add(cell, 2.0 * gain),
        );
        for got in removed_both_ways(&model, &block, &edges, r) {
            assert!(
                (got - want).abs() < 3e-6,
                "draft {draft}: {got} against {want}"
            );
        }
    }
}

/// The runs two lists of disjoint runs share.
fn meet(a: &[(f64, f64)], b: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for &(p, q) in a {
        for &(r, t) in b {
            let (lo, hi) = (p.max(r), q.min(t));
            if hi > lo {
                out.push((lo, hi));
            }
        }
    }
    out
}

/// A block eight square and two high notched at the top corner `(4, 3, 2)`
/// past two walls: one through the crease `y = 3` (the first edge, along
/// `x` from the corner), and one along the second edge, which leaves the
/// corner at `plan` from the first in plan. Each wall leans from the
/// vertical by its own angle, its foot in under the material for a
/// positive lean. Returns the block and the two top creases leaving the
/// corner, the one along `x` first.
fn turned_corner(
    model: &mut ogeom_topo::Model,
    plan: f64,
    [lean, leant]: [f64; 2],
) -> (ogeom_topo::Shape, [ogeom_topo::Shape; 2]) {
    let block = ogeom_algo::make_box(model, Frame::WORLD, (8.0, 8.0, 2.0), T)
        .unwrap()
        .shape;
    let corner = Point::new(4.0, 3.0, 2.0);
    let (s, c) = plan.sin_cos();
    let past = |model: &mut ogeom_topo::Model, n: Vector, x: Vector| {
        let y = n.cross(x);
        let origin = corner - x * 20.0 - y * 20.0;
        let frame = Frame::new(
            origin,
            Direction::new(n, T).unwrap(),
            Direction::new(x, T).unwrap(),
            T,
        )
        .unwrap();
        ogeom_algo::make_box(model, frame, (40.0, 40.0, 40.0), T)
            .unwrap()
            .shape
    };
    let first = past(model, Vector::new(0.0, lean.cos(), -lean.sin()), Vector::X);
    let second = past(
        model,
        Vector::new(s * leant.cos(), -c * leant.cos(), -leant.sin()),
        Vector::new(c, s, 0.0),
    );
    let notch = ogeom_bool::common(model, &first, &second, T).unwrap().shape;
    let cut = ogeom_bool::cut(model, &block, &notch, T).unwrap().shape;
    let block = explore_unique(model, &cut, ShapeType::Solid).unwrap()[0].clone();
    let point = |v: &ogeom_topo::Shape| {
        let p = model.node(v).unwrap().data().as_vertex().unwrap().point;
        v.transform(model.datums()).unwrap().apply(p)
    };
    // The top creases leaving the corner, the one along `x` first.
    let mut edges: Vec<ogeom_topo::Shape> = explore_unique(model, &block, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(model, e).unwrap().unwrap();
            let (p, q) = (point(&a), point(&b));
            (p.z - 2.0).abs() < 1e-9
                && (q.z - 2.0).abs() < 1e-9
                && (p.distance(corner) < 1e-9 || q.distance(corner) < 1e-9)
        })
        .collect();
    edges.sort_by_key(|e| {
        let (a, b) = ogeom_algo::edge_vertices(model, e).unwrap().unwrap();
        (point(&a).y - 3.0).abs() + (point(&b).y - 3.0).abs() > 1e-9
    });
    let [first, second] = edges.as_slice() else {
        panic!("two top creases at the corner");
    };
    (block, [first.clone(), second.clone()])
}

/// Two top creases rounded where they meet at an inside corner off square
/// in plan, or with both walls leaning by different angles. The walls'
/// meeting line leans, so material beside each band past the other wall's
/// plane, before that band's cap, lies above both bands; and the two
/// bands' widenings meet past both walls, where the corner is open. Each
/// order takes what lies above both bands carried on, measured within a
/// box round the corner a column at a time. Capped at the vertex instead,
/// the bands left the corner standing.
#[test]
fn rounds_at_inside_corners_off_square_or_leaning_both_ways_take_all_above_both() {
    let r = 0.5;
    let lean = 0.1_f64.atan();
    let (wx, wy, half) = (4.0, 3.0, 1.2);
    for (degrees, leans) in [
        (120.0_f64, [lean, 0.0]),
        (150.0, [lean, 0.0]),
        (170.0, [lean, 0.0]),
        (90.0, [lean, 3.0_f64.to_radians()]),
        (90.0, [-lean, (-2.0_f64).to_radians()]),
        (120.0, [lean, 3.0_f64.to_radians()]),
        (150.0, [-lean, 2.0_f64.to_radians()]),
    ] {
        let plan = degrees.to_radians();
        let (s, c) = plan.sin_cos();
        let first = Section::new(0.0, (leans[0].cos(), -leans[0].sin()), r);
        let second = Section::new(0.0, (leans[1].cos(), -leans[1].sin()), r);
        let want = integral(
            (wx - half, wx + half),
            (wy - half, wy + half),
            2e3,
            |x, y| {
                let (px, py) = (x - wx, y - wy);
                let (ua, ub) = (py, px * s - py * c);
                let [a0, a1] = first.column(ua);
                let [b0, b1] = second.column(ub);
                let (pa, pb) = (first.past(ua, 0.0), second.past(ub, 0.0));
                // Under the top, out of the open corner past both walls.
                let above = meet(&[a0, a1, pa], &[b0, b1, pb]);
                let runs: Vec<(f64, f64)> = above
                    .iter()
                    .map(|&(p, q)| (p.max(0.0), q.min(2.0)))
                    .collect();
                span(&runs) - overlap(&runs, &meet(&[pa], &[pb]))
            },
        );
        let mut model = ogeom_topo::Model::new();
        let (block, edges) = turned_corner(&mut model, plan, leans);
        for named in [edges.clone(), [edges[1].clone(), edges[0].clone()]] {
            let mut model = model.clone();
            let rounded = ogeom_fillet::fillet_edges(&mut model, &block, &named, r, T).unwrap();
            assert!(
                ogeom_algo::check(&model, &rounded.shape, T)
                    .unwrap()
                    .is_valid()
            );
            let gone = ogeom_bool::cut(&mut model, &block, &rounded.shape, T)
                .unwrap()
                .shape;
            let window = ogeom_algo::make_box(
                &mut model,
                Frame::new(
                    Point::new(wx - half, wy - half, -1.0),
                    Direction::Z,
                    Direction::X,
                    T,
                )
                .unwrap(),
                (half * 2.0, half * 2.0, 4.0),
                T,
            )
            .unwrap()
            .shape;
            let near = ogeom_bool::common(&mut model, &gone, &window, T)
                .unwrap()
                .shape;
            let got = volume(&model, &near, 1e-4);
            assert!(
                (got - want).abs() < 5e-6,
                "{degrees} degrees, leaning {leans:?}: {got} against {want}"
            );
        }
    }
}

/// A box's top edges rounded, two adjacent ones and all four: each band
/// meets its neighbours in the plane bisecting their square outside
/// corner, and each band's face stays credited to its own edge. Each takes
/// the corner outside its quarter circle along the box less, at each
/// corner, the cell past the bisecting plane its cap overlapped the other
/// band by, the corner's moment about its wall: (1 - pi/4) r^2 L per edge
/// less 2 r^3 (5/6 - pi/4) per corner.
#[test]
fn a_box_s_top_edges_meet_at_their_corners_each_credited_its_own_band() {
    let pi = core::f64::consts::PI;
    let r = 0.5;
    let section = (1.0 - pi / 4.0) * r * r;
    let cell = 2.0 * (5.0 / 6.0 - pi / 4.0) * r.powi(3);
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (4.0, 4.0, 2.0), T)
        .unwrap()
        .shape;
    let top = [(2.0, 0.0), (4.0, 2.0), (2.0, 4.0), (0.0, 2.0)]
        .map(|(x, y)| placed_edge_through(&model, &block, Point::new(x, y, 2.0)));
    for (count, corners) in [(2_usize, 1.0), (4, 4.0)] {
        let named = &top[..count];
        let removed = removed_both_ways(&model, &block, named, r);
        #[allow(clippy::cast_precision_loss, reason = "a handful of edges")]
        let want = section.mul_add(4.0 * count as f64, -cell * corners);
        for got in removed {
            assert!(
                (got - want).abs() < 2e-5,
                "{count} edges: {got} against {want}"
            );
        }
        let mut rounded_model = model.clone();
        let rounded = ogeom_fillet::fillet_edges(&mut rounded_model, &block, named, r, T).unwrap();
        for edge in named {
            assert_eq!(
                rounded.history.generated(edge).len(),
                1,
                "{count} edges: each edge's own band"
            );
        }
    }
}

/// A prism two high whose top's outline runs from `(8, 3)` to the corner
/// `(4, 3)` along an arc bulging 0.4 towards +y, then turns `turn` back
/// from the arc's end tangent along a straight edge three long. Returns
/// the solid and the arc's and the straight edge's top creases.
fn curved_turn(
    model: &mut ogeom_topo::Model,
    turn: f64,
) -> (ogeom_topo::Shape, [ogeom_topo::Shape; 2]) {
    use ogeom_geom::{Curve, Curve3d as _, LineCurve};
    use ogeom_math::Axis;
    let (from, corner) = (Point::new(8.0, 3.0, 0.0), Point::new(4.0, 3.0, 0.0));
    let (half, sag): (f64, f64) = (2.0, 0.4);
    let radius = half.mul_add(half, sag * sag) / (2.0 * sag);
    let centre = Point::new(6.0, 3.0 + sag - radius, 0.0);
    let frame = Frame::new(
        centre,
        Direction::Z,
        Direction::new(from - centre, T).unwrap(),
        T,
    )
    .unwrap();
    let arc: Curve = CircleCurve::new(Circle::new(frame, radius, T).unwrap()).into();
    let sweep = {
        let (u, v) = ((from - centre) / radius, (corner - centre) / radius);
        u.cross(v).z.atan2(u.dot(v))
    };
    let ahead = {
        let d = arc.d1_at(sweep, T).unwrap();
        d / d.magnitude()
    };
    let (s, c) = turn.sin_cos();
    let along = Vector::new(
        ahead.x.mul_add(c, ahead.y * s),
        (-ahead.x).mul_add(s, ahead.y * c),
        0.0,
    );
    let tip = corner + along * 3.0;
    let outline = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(8.0, 0.0, 0.0),
        from,
        corner,
        tip,
        Point::new(tip.x, 9.0, 0.0),
        Point::new(0.0, 9.0, 0.0),
    ];
    let vertices: Vec<ogeom_topo::Shape> = outline
        .iter()
        .map(|p| ogeom_algo::make_vertex(model, *p).shape)
        .collect();
    let mut edges = Vec::new();
    for i in 0..outline.len() {
        let j = (i + 1) % outline.len();
        let (p, q) = (outline[i], outline[j]);
        let edge = if i == 2 {
            ogeom_algo::make_edge_between(
                model,
                arc.clone(),
                (0.0, sweep),
                &vertices[i],
                &vertices[j],
                T,
            )
        } else {
            let line: Curve =
                LineCurve::new(Axis::new(p, Direction::new(q - p, T).unwrap())).into();
            ogeom_algo::make_edge_between(
                model,
                line,
                (0.0, p.distance(q)),
                &vertices[i],
                &vertices[j],
                T,
            )
        };
        edges.push(edge.unwrap().shape);
    }
    let wire = ogeom_algo::make_wire(model, &edges, T).unwrap().shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let prism = ogeom_algo::make_prism(model, &face, Vector::new(0.0, 0.0, 2.0), T)
        .unwrap()
        .shape;
    let lift = Vector::new(0.0, 0.0, 2.0);
    let mid_arc = arc.point_at(sweep / 2.0, T).unwrap() + lift;
    let along_arc = explore_unique(model, &prism, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let b = ogeom_algo::shape_bounds(model, e, T).unwrap();
            let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
            (lo.z - 2.0).abs() < 1e-6
                && (hi.z - 2.0).abs() < 1e-6
                && hi.y >= mid_arc.y - 1e-6
                && (lo.x - 4.0).abs() < 1e-6
                && (hi.x - 8.0).abs() < 1e-6
        })
        .expect("the arc's top crease");
    let straight = placed_edge_through(model, &prism, corner + along * 1.5 + lift);
    (prism, [along_arc, straight])
}

/// Two top creases, an arc's and a straight edge's, rounded where the
/// outline turns back on itself past them. At 150 degrees the corner cell
/// both rounds hold lies in the material and the two meet in the plane
/// bisecting their straight continuations, in either order. Turned 170 or
/// 178 degrees back, the continuations reach that plane only past the
/// solid's far side: no junction is taken, and the request rounds exactly
/// as it does with none. Those figures are what the rounding without
/// junctions takes.
#[test]
fn rounds_at_a_curved_crease_turning_nearly_back_on_itself() {
    let mut model = ogeom_topo::Model::new();
    let (prism, edges) = curved_turn(&mut model, 150.0_f64.to_radians());
    let removed = removed_both_ways(&model, &prism, &edges, 0.5);
    assert!((removed[0] - removed[1]).abs() < 2e-6, "{removed:?}");
    assert!((removed[0] - 0.429_067).abs() < 2e-6, "{removed:?}");
    for (degrees, want) in [
        (170.0_f64, [0.379_278, 0.380_335]),
        (178.0, [0.379_036, 0.379_246]),
    ] {
        let mut model = ogeom_topo::Model::new();
        let (prism, edges) = curved_turn(&mut model, degrees.to_radians());
        let removed = removed_both_ways(&model, &prism, &edges, 0.5);
        for (got, want) in removed.iter().zip(want) {
            assert!(
                (got - want).abs() < 2e-6,
                "{degrees} degrees: {removed:?} against {want}"
            );
        }
    }
}

/// A block eight square and two high with a pocket `[2, 6]` square and one
/// deep whose wall `y = 2` leans `drafts[0]` per unit of depth into the
/// pocket's floor and whose wall `x = 6` leans `drafts[1]`, all turned by
/// `turn`: the block and the pocket's four rim creases, along `y = 2`
/// first, then round the rim.
fn leaning_pocket(
    model: &mut ogeom_topo::Model,
    drafts: [f64; 2],
    turn: &ogeom_math::Transform,
) -> (ogeom_topo::Shape, Vec<ogeom_topo::Shape>) {
    let block = ogeom_algo::make_box(model, Frame::WORLD, (8.0, 8.0, 2.0), T)
        .unwrap()
        .shape;
    // The half-space past the plane through `on` with unit normal `n`.
    let past = |model: &mut ogeom_topo::Model, on: Point, n: Vector| {
        let seed = if n.x.abs() < 0.9 {
            Vector::X
        } else {
            Vector::Y
        };
        let x = seed - n * seed.dot(n);
        let x = x / x.magnitude();
        let y = n.cross(x);
        let frame = Frame::new(
            on - x * 20.0 - y * 20.0,
            Direction::new(n, T).unwrap(),
            Direction::new(x, T).unwrap(),
            T,
        )
        .unwrap();
        ogeom_algo::make_box(model, frame, (40.0, 40.0, 40.0), T)
            .unwrap()
            .shape
    };
    let unit = |v: Vector| v / v.magnitude();
    let top = |x: f64, y: f64| Point::new(x, y, 2.0);
    let sides = [
        past(model, top(4.0, 2.0), unit(Vector::new(0.0, 1.0, drafts[0]))),
        past(
            model,
            top(6.0, 4.0),
            unit(Vector::new(-1.0, 0.0, drafts[1])),
        ),
        past(model, top(4.0, 6.0), Vector::new(0.0, -1.0, 0.0)),
        past(model, top(2.0, 4.0), Vector::X),
        past(model, Point::new(4.0, 4.0, 1.0), Vector::Z),
    ];
    let mut pocket = sides[0].clone();
    for side in &sides[1..] {
        pocket = ogeom_bool::common(model, &pocket, side, T).unwrap().shape;
    }
    let cut = ogeom_bool::cut(model, &block, &pocket, T).unwrap().shape;
    let block = explore_unique(model, &cut, ShapeType::Solid).unwrap()[0].clone();
    let block = model.placed(&block, *turn);
    let edges = [(4.0, 2.0), (6.0, 4.0), (4.0, 6.0), (2.0, 4.0)]
        .map(|(x, y)| placed_edge_through(model, &block, turn.apply(top(x, y))))
        .to_vec();
    (block, edges)
}

/// A square pocket's four rim creases rounded together: each inside corner
/// meets flush in its bisecting plane, the junction closing the loop as
/// well as the three before it, so the four bands join into one ring.
/// Each crease takes its section along its four, and each corner the cell
/// of its two quarter rounds, 2 r^3 (5/6 - pi/4), in any order the edges
/// are named in. The junction closing the loop was dropped, and two of
/// the corners kept their cells.
#[test]
fn a_square_pocket_s_rim_meets_at_every_corner_in_any_order() {
    let r: f64 = 0.5;
    let pi = core::f64::consts::PI;
    let want = 4.0f64.mul_add(
        4.0 * (1.0 - pi / 4.0) * r * r,
        8.0 * (5.0 / 6.0 - pi / 4.0) * r.powi(3),
    );
    let mut model = ogeom_topo::Model::new();
    let (block, edges) = drafted_pocket(&mut model, 0.0);
    for order in [[0, 1, 2, 3], [2, 0, 3, 1], [1, 3, 0, 2]] {
        let named: Vec<ogeom_topo::Shape> = order.iter().map(|&k| edges[k].clone()).collect();
        for got in removed_both_ways(&model, &block, &named, r) {
            assert!((got - want).abs() < 2e-6, "{order:?}: {got} against {want}");
        }
    }
}

/// A pocket whose two adjacent walls lean by different drafts, rounded
/// round its rim, as built and turned 0.7 about (1, 1, 1): its two drafted
/// corners and the one between the drafted walls are crossed, and the part
/// turned takes what it takes as built, in either order. Turned, a band's
/// end arc came out of the boolean in two pieces at its circle's seam,
/// which follows the placement, and the corner was refused.
#[test]
fn a_leaning_pocket_s_rim_rounds_alike_however_turned() {
    let spin = ogeom_math::Transform::rotation(
        ogeom_math::Axis::new(
            Point::ORIGIN,
            Direction::new(Vector::new(1.0, 1.0, 1.0), T).unwrap(),
        ),
        0.7,
    );
    let mut removed = Vec::new();
    for turn in [ogeom_math::Transform::IDENTITY, spin] {
        let mut model = ogeom_topo::Model::new();
        let (block, edges) = leaning_pocket(&mut model, [0.1, 0.05], &turn);
        removed.extend(removed_both_ways(&model, &block, &edges, 0.5));
    }
    for got in &removed {
        assert!((got - removed[0]).abs() < 1e-6, "{removed:?}");
    }
    // Taken crossed, past what the corners capped leave.
    assert!(removed[0] > 0.83, "{removed:?}");
}

/// What a square inside corner between two walls leaning by `leans` (as
/// [`turned_corner`] leans them) takes beyond its two creases' sections
/// along their lengths: what lies above both bands carried on within a
/// box round the corner, a column at a time, less the two sections along
/// the box's reach.
fn leaning_corner_gain(leans: [f64; 2], r: f64) -> f64 {
    let half = 1.2;
    let first = Section::new(0.0, (leans[0].cos(), -leans[0].sin()), r);
    let second = Section::new(0.0, (leans[1].cos(), -leans[1].sin()), r);
    let window = integral((-half, half), (-half, half), 2e3, |x, y| {
        let [a0, a1] = first.column(y);
        let [b0, b1] = second.column(x);
        let (pa, pb) = (first.past(y, 0.0), second.past(x, 0.0));
        let above = meet(&[a0, a1, pa], &[b0, b1, pb]);
        let runs: Vec<(f64, f64)> = above
            .iter()
            .map(|&(p, q)| (p.max(0.0), q.min(2.0)))
            .collect();
        span(&runs) - overlap(&runs, &meet(&[pa], &[pb]))
    });
    half.mul_add(-(first.area + second.area), window)
}

/// A block eight square and two high with a pocket over `[x0, x1]` by
/// `[y0, y1]`, one deep, each wall leaning `d` per unit of depth in turn
/// round from `y = y0` (a positive draft narrowing the pocket's floor):
/// the block and the pocket's four rim creases in that order.
fn pocket_of(
    model: &mut ogeom_topo::Model,
    [x0, y0, x1, y1]: [f64; 4],
    d: [f64; 4],
) -> (ogeom_topo::Shape, Vec<ogeom_topo::Shape>) {
    let block = ogeom_algo::make_box(model, Frame::WORLD, (8.0, 8.0, 2.0), T)
        .unwrap()
        .shape;
    let past = |model: &mut ogeom_topo::Model, on: Point, n: Vector| {
        let n = n / n.magnitude();
        let seed = if n.x.abs() < 0.9 {
            Vector::X
        } else {
            Vector::Y
        };
        let x = seed - n * seed.dot(n);
        let x = x / x.magnitude();
        let y = n.cross(x);
        let frame = Frame::new(
            on - x * 20.0 - y * 20.0,
            Direction::new(n, T).unwrap(),
            Direction::new(x, T).unwrap(),
            T,
        )
        .unwrap();
        ogeom_algo::make_box(model, frame, (40.0, 40.0, 40.0), T)
            .unwrap()
            .shape
    };
    let (xm, ym) = ((x0 + x1) / 2.0, (y0 + y1) / 2.0);
    let top = |x: f64, y: f64| Point::new(x, y, 2.0);
    let sides = [
        past(model, top(xm, y0), Vector::new(0.0, 1.0, d[0])),
        past(model, top(x1, ym), Vector::new(-1.0, 0.0, d[1])),
        past(model, top(xm, y1), Vector::new(0.0, -1.0, d[2])),
        past(model, top(x0, ym), Vector::new(1.0, 0.0, d[3])),
        past(model, Point::new(xm, ym, 1.0), Vector::Z),
    ];
    let mut pocket = sides[0].clone();
    for s in &sides[1..] {
        pocket = ogeom_bool::common(model, &pocket, s, T).unwrap().shape;
    }
    let cut = ogeom_bool::cut(model, &block, &pocket, T).unwrap().shape;
    let block = explore_unique(model, &cut, ShapeType::Solid).unwrap()[0].clone();
    let edges = [(xm, y0), (x1, ym), (xm, y1), (x0, ym)]
        .map(|(x, y)| placed_edge_through(model, &block, top(x, y)))
        .to_vec();
    (block, edges)
}

/// Pockets whose walls lean by different drafts rounded round their
/// rims: four walls each its own way, and a pocket a millimetre and a half
/// wide whose two long walls lean. Every corner is crossed or flush, each
/// band crossed at both ends cut in two, and each crease takes its section
/// along its length and each corner what lies above both its bands beyond
/// those sections, in either order.
#[test]
fn pockets_whose_walls_lean_their_own_ways_meet_at_every_corner() {
    let r = 0.5;
    let lean = |d: f64| -d.atan();
    let area = |d: f64| Section::new(0.0, (lean(d).cos(), -lean(d).sin()), r).area;
    for (rect, d) in [
        ([2.0, 2.0, 6.0, 6.0], [0.1, 0.05, -0.05, 0.08]),
        ([2.0, 3.25, 6.0, 4.75], [0.1, 0.0, 0.05, 0.0]),
    ] {
        let mut model = ogeom_topo::Model::new();
        let (block, edges) = pocket_of(&mut model, rect, d);
        let (w, h) = (rect[2] - rect[0], rect[3] - rect[1]);
        let mut want = w.mul_add(area(d[0]) + area(d[2]), h * (area(d[1]) + area(d[3])));
        for (a, b) in [(0, 3), (0, 1), (2, 1), (2, 3)] {
            want += leaning_corner_gain([lean(d[a]), lean(d[b])], r);
        }
        // Each corner's box is integrated to about a micron cubed.
        for got in removed_both_ways(&model, &block, &edges, r) {
            assert!((got - want).abs() < 1e-5, "{d:?}: {got} against {want}");
        }
    }
}
