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
/// past the corner's vertex, within the vertex's tolerance. Each blend's
/// cap stands at the vertex, where the other blend is flush with it, so
/// each removes exactly its edge's section over the edge's length from the
/// outside wall to the inside corner, in either order. A cap left where
/// the curve ends would cut into the other blend past the corner.
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
            let want = 1500.0 - section * 20.0;
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
/// that wall, the first does. Each blend's cap stands where its line meets
/// the other edge's, flush with the other blend, so each removes exactly
/// its edge's section over the edge's length, in either order. A cap at
/// the vertex would cut into the other blend past the corner.
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
            let want = 1500.0 - section * 20.0;
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
