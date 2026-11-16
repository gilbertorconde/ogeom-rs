//! Filleting a tangent chain of edges in one call.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/swept.rs"]
#[allow(
    dead_code,
    reason = "this file sums its oracles with the adaptive rule alone"
)]
mod swept;

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

/// Every edge of an L-bracket, in either order. The step's inner corners,
/// where its re-entrant edge meets two convex ones, are no convex corners:
/// each closes with one setback patch, and the two orders give the same
/// solid.
#[test]
fn every_edge_of_an_l_bracket_rounds_in_either_order() {
    let mut volumes = Vec::new();
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
        // Nothing of it outside the bracket and the re-entrant edge's
        // concave fill: the convex rounds only take material off.
        let deflection = ogeom_mesh::Deflection {
            chord: 1e-4,
            ..ogeom_mesh::Deflection::default()
        };
        let mut past = f64::NEG_INFINITY;
        for f in explore_unique(&model, &rounded, ShapeType::Face).unwrap() {
            for p in ogeom_mesh::triangulate_face(&model, &f, deflection, T)
                .unwrap()
                .positions
            {
                let in_box = [-p.x, -p.y, -p.z, p.x - 2.0, p.y - 2.0, p.z - 2.0]
                    .into_iter()
                    .fold(f64::NEG_INFINITY, f64::max);
                let bracket = in_box.max((p.x - 1.0).min(p.z - 1.0));
                let fill = [1.0 - p.x, p.x - 1.1, 1.0 - p.z, p.z - 1.1, -p.y, p.y - 2.0]
                    .into_iter()
                    .fold(0.1 - (p.x - 1.1).hypot(p.z - 1.1), f64::max);
                past = past.max(bracket.min(fill));
            }
        }
        assert!(
            past <= 2e-4,
            "reversed {reversed}: {past} outside the bracket and its fill"
        );
        let v = volume(&model, &rounded, 1e-4);
        assert!((v - 6.0).abs() < 0.1, "reversed {reversed}: {v}");
        volumes.push(v.to_bits());
    }
    assert_eq!(
        volumes[0], volumes[1],
        "the two orders round different solids"
    );
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
/// degrees. The rim stops where the neck leaves it, at an inside corner
/// whose third edge, the neck's top crease, is left sharp: there the band
/// ends flush on the neck's wall (ADR 0001, amendment "M5e", a free end;
/// SMLib's mixed 3x1 corner, "the end face will be extended out to the
/// fillet surface"), so near that end the part of the band past the
/// wall's plane stays. Over the seam it ends in the plane through the axis
/// (the long arc unasked) or runs on into the long arc's band. Rounded
/// alone, and with the long arc, each round removes the corner between
/// the top and the wall turned about the axis, a column at radius `rho`
/// taken over the turn where it stands on the crease's side of the neck's
/// walls `y = -1` and `y = 1` ([`keyhole_removed`]). The current path
/// capped the band in the plane through the axis at the wall, the whole
/// corner over 20 and 300 degrees by Pappus; that figure was its own.
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
    let before = volume(&model, &keyhole, 1e-4);
    // The turn a column at `rho` is taken over: from the seam to where it
    // crosses the wall `y = -1`, and with the long arc from where it
    // crosses the wall `y = 1`.
    let alone = |rho: f64| 2.0 * pi - (1.0 / rho).asin() - 310_f64.to_radians();
    let both = |rho: f64| 2.0 * pi - 2.0 * (1.0 / rho).asin();
    for (named, sweep) in [
        (vec![short.clone()], keyhole_removed(radius, r, alone)),
        (vec![long, short], keyhole_removed(radius, r, both)),
    ] {
        let mut model = model.clone();
        let rounded = ogeom_fillet::fillet_edges(&mut model, &keyhole, &named, r, T).unwrap();
        assert!(
            ogeom_algo::check(&model, &rounded.shape, T)
                .unwrap()
                .is_valid()
        );
        let removed = before - volume(&model, &rounded.shape, 1e-4);
        let want = sweep;
        assert!(
            (removed - want).abs() < want * 1e-4,
            "{} arcs: {removed} against {want}",
            named.len()
        );
    }
}

/// What a rim's round of radius `r` on a drum of radius `big` takes off a
/// top `z = 1` over the turn `turn(rho)` a column at radius `rho` is taken
/// over: the column's height under the top, `r - sqrt(r^2 - (rho - big +
/// r)^2)` outside the ball, times `rho` times the turn, integrated across
/// the band, `rho = big - r + r sin phi` taking the square root's end out.
fn keyhole_removed(big: f64, r: f64, turn: impl Fn(f64) -> f64) -> f64 {
    gauss(
        |phi| {
            let rho = r.mul_add(phi.sin(), big - r);
            let height = r * (1.0 - phi.cos());
            height * rho * turn(rho) * r * phi.cos()
        },
        0.0,
        core::f64::consts::FRAC_PI_2,
    )
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
/// pi/4). The bevels mitre across the corner the same way: the top's
/// outline set back round its reflex corner, each bevel running on past the
/// vertex to the other's plane, which takes the column of the corner square
/// both strips reach, `d³/3` (the integral of `s²` over the setback `s`
/// from 0 to `d`). A cap left where the curve ends would cut into the
/// other blend past the corner.
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
                size * size * size / 3.0
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
/// the bevels running on past it to each other's plane (the column `d³/3`,
/// as above, read on the faces' exact planes whatever the vertex does), and
/// the rounds in the plane bisecting the corner through it, as above, in
/// either order. A cap at the vertex would cut into the other blend past
/// the corner.
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
            // The faces are exact planes: the bevels are the plain L's,
            // each running on past the vertex to the other's plane at the
            // reflex corner and taking the column `d³/3` both strips reach,
            // as above. (The wedges capped them at the lines' meeting and
            // left that column standing.)
            let cell = if chamfer {
                size * size * size / 3.0
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

/// Both orders of rounding the two top creases of [`turned_corner`] at
/// radius `r`, each checked valid and measured within a box round the
/// corner, against what lies above both bands carried on there, a column
/// at a time.
fn takes_all_above_both(degrees: f64, leans: [f64; 2], r: f64) {
    let (wx, wy, half) = (4.0, 3.0, 1.2);
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
    let lean = 0.1_f64.atan();
    for (degrees, leans) in [
        (120.0_f64, [lean, 0.0]),
        (150.0, [lean, 0.0]),
        (170.0, [lean, 0.0]),
        (90.0, [lean, 3.0_f64.to_radians()]),
        (90.0, [-lean, (-2.0_f64).to_radians()]),
        (120.0, [lean, 3.0_f64.to_radians()]),
        (150.0, [-lean, 2.0_f64.to_radians()]),
    ] {
        takes_all_above_both(degrees, leans, 0.5);
    }
}

/// Two top creases rounded where they meet at an inside corner whose walls
/// lean alike, their feet out under the corner as a moulded pocket's are:
/// the two end sections are reflections of each other in the plane
/// bisecting the turn, so the bands meet flush there, but each band's end
/// section reaches past that plane down its wall. Each band is the
/// straight band of its end section on its own side of the plane, before
/// its cap as past it, and the two take what lies above both, square in
/// plan, off square and acute. Carried on from its cap alone, the band
/// could not reach back to the plane, and both were capped at the vertex.
#[test]
fn rounds_at_inside_corners_whose_walls_lean_alike_under_them_meet_in_the_bisecting_plane() {
    let lean = -0.1_f64.atan();
    for (degrees, leans) in [
        (90.0_f64, [lean, lean]),
        (120.0, [lean, lean]),
        (75.0, [lean, lean]),
        (60.0, [-0.05_f64.atan(), -0.05_f64.atan()]),
    ] {
        takes_all_above_both(degrees, leans, 0.5);
    }
}

/// A box's top edges rounded, two adjacent ones and all four: each band
/// meets its neighbours in the plane bisecting their square outside
/// corner, and each band's face stays credited to its own edge. Each takes
/// the corner outside its quarter circle along the box less, at each
/// corner, the cell past the bisecting plane its cap overlapped the other
/// band by, the corner's moment about its wall: (1 - pi/4) r^2 L per edge
/// less 2 r^3 (5/6 - pi/4) per corner. The solid keeps no vertex but the
/// box's own, the bands' ends and the mitres': the circle of a band's end
/// arc crosses the next band's rail again a diameter past the corner,
/// outside every face, and splits nothing there.
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
    // Two adjacent edges keep the box's four bottom corners, its one top
    // corner off both bands, two band ends at each free end and the rail
    // ends and mitre's top at the shared corner, and Euler's 2 gives the
    // edges; all four keep the bottom corners, the rail ends and the
    // mitres' tops.
    for (count, corners, minimal) in [(2_usize, 1.0, (8, 17, 11)), (4, 4.0, (10, 20, 12))] {
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
        let diagnosis = ogeom_algo::check(&rounded_model, &rounded.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{count} edges: {diagnosis}");
        let held = |kind| {
            explore_unique(&rounded_model, &rounded.shape, kind)
                .unwrap()
                .len()
        };
        assert_eq!(
            (
                held(ShapeType::Face),
                held(ShapeType::Edge),
                held(ShapeType::Vertex)
            ),
            minimal,
            "{count} edges: faces, edges and vertices"
        );
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

/// How a notch's corner is closed in [`turn_removed_by`]: the bands
/// carried on through the inside corner (amendment "M5e"), the ball rolled
/// round the notch's edge, the bands alone with nothing between them, or
/// the bands with the whole sector behind the corner taken to the depth
/// where the ball touches the edge.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Corner {
    Joined,
    Rolled,
    Flat,
    Hollow,
}

/// What rounding [`curved_turn`]'s two top creases at `r` takes, summed
/// independently of any boolean, the corner closed as `corner` says.
///
/// [`Corner::Joined`], by the rule of the bands carried on through an
/// inside corner (ADR 0001, amendment "M5e"): the arc's band carried on is
/// its torus band over the arc, then each end's section run straight on
/// along the arc's tangent there; the straight edge's is its band along the
/// whole line. Past the arc's wall is outside its cylinder, past the
/// straight wall the open side of its plane. A point of the solid goes
/// where it lies in one band and past the other's wall or in the other's
/// band.
///
/// [`Corner::Rolled`], by the rule for a notch whose rails meet only past
/// the top (ADR 0001, amendment "corners where a round's rail misses"):
/// each band ends where its ball touches the notch's edge, in the plane
/// through that edge square to its own crease, and the ball rolls round the
/// edge behind the corner, touching it and the top. At each height that
/// roll takes the sector of the disc of the band's width about the corner,
/// through the turn, centred on the way back from the notch. [`Corner::Flat`]
/// and [`Corner::Hollow`] bound it: nothing taken there, or the sector of
/// the radius at every height down to the ball's touch on the edge.
///
/// At height `z` a band is the strip `w = r - sqrt(r^2 - t^2)` in from its
/// wall, `t = z - (2 - r)`, so each height's slice is integrated across
/// `x`, each column's length in `y` read exactly between the roots of every
/// boundary there, and the slices summed over `t = r sin(phi)`.
fn turn_removed_by(turn: f64, r: f64, corner_by: Corner) -> f64 {
    turn_removed_within(turn, r, corner_by, false)
}

/// [`turn_removed_by`], with `sector`: only what lies in the sector of the
/// radius behind the notch's corner, through the turn, centred on the way
/// back from the notch (the region the rolled corner stands in).
fn turn_removed_within(turn: f64, r: f64, corner_by: Corner, sector: bool) -> f64 {
    let pivot = corner_by != Corner::Joined;
    let (cx, cy, big) = (6.0, 3.4 - 5.2, 5.2);
    let angle = |x: f64, y: f64| (y - cy).atan2(x - cx);
    let (tf, tc) = (angle(8.0, 3.0), angle(4.0, 3.0));
    let at = |t: f64| (cx + big * t.cos(), cy + big * t.sin());
    // Each end, its tangent running on past it and its outward radial.
    let ends = [
        (at(tf), (tf.sin(), -tf.cos()), (tf.cos(), tf.sin())),
        (at(tc), (-tc.sin(), tc.cos()), (tc.cos(), tc.sin())),
    ];
    let ahead = ends[1].1;
    let (s, c) = turn.sin_cos();
    let along = (
        ahead.0.mul_add(c, ahead.1 * s),
        (-ahead.0).mul_add(s, ahead.1 * c),
    );
    let into = (-along.1, along.0);
    let corner = (4.0, 3.0);
    let tip = (corner.0 + 3.0 * along.0, corner.1 + 3.0 * along.1);
    let outline = [
        (0.0, 0.0),
        (8.0, 0.0),
        (8.0, 3.0),
        corner,
        tip,
        (tip.0, 9.0),
        (0.0, 9.0),
    ];
    let dot = |a: (f64, f64), b: (f64, f64)| a.0.mul_add(b.0, a.1 * b.1);
    let rho = |x: f64, y: f64| (x - cx).hypot(y - cy);
    let in_outline = |x: f64, y: f64| {
        let mut inside = false;
        for k in 0..outline.len() {
            let (a, b) = (outline[k], outline[(k + 1) % outline.len()]);
            if (a.1 > y) != (b.1 > y) && x < a.0 + (b.0 - a.0) * (y - a.1) / (b.1 - a.1) {
                inside = !inside;
            }
        }
        // The chord's polygon, or the segment between the chord and the arc.
        inside || (y > 3.0 && rho(x, y) < big)
    };
    // The way back from the notch, and the cosine of half the turn.
    let back = {
        let b = (-(-ahead.0 + along.0), -(-ahead.1 + along.1));
        let l = b.0.hypot(b.1);
        (b.0 / l, b.1 / l)
    };
    let half = (0.5 * turn).cos();
    let in_sector = |x: f64, y: f64| {
        let off = (x - corner.0, y - corner.1);
        let reach = off.0.hypot(off.1);
        reach <= r && (reach == 0.0 || dot(off, back) >= half * reach)
    };
    let removed_at = |w: f64, x: f64, y: f64| {
        if !in_outline(x, y) || (sector && !in_sector(x, y)) {
            return false;
        }
        if pivot {
            let off = (x - corner.0, y - corner.1);
            let d = dot(off, into);
            let b_band = (0.0..=w).contains(&d) && dot(off, along) >= 0.0;
            let q = rho(x, y);
            let torus = (big - w..=big).contains(&q) && (tf..=tc).contains(&angle(x, y));
            let (e, t, n) = ends[0];
            let far = (x - e.0, y - e.1);
            let run_on = dot(far, t) >= 0.0 && (0.0..=w).contains(&-dot(far, n));
            let reach = off.0.hypot(off.1);
            let depth = match corner_by {
                Corner::Flat => -1.0,
                Corner::Hollow => r,
                _ => w,
            };
            let rolled = reach <= depth && (reach == 0.0 || dot(off, back) >= half * reach);
            return torus || run_on || b_band || rolled;
        }
        let d = dot((x - corner.0, y - corner.1), into);
        let in_b = (0.0..=w).contains(&d);
        let past_b = d < 0.0;
        let q = rho(x, y);
        let torus = (big - w..=big).contains(&q) && (tf..=tc).contains(&angle(x, y));
        let run_on = ends.iter().any(|(e, t, n)| {
            let off = (x - e.0, y - e.1);
            dot(off, *t) >= 0.0 && (0.0..=w).contains(&-dot(off, *n))
        });
        let in_a = torus || run_on;
        let past_a = q > big;
        (in_a && (past_b || in_b)) || (in_b && (past_a || in_a))
    };
    // A column's length in `y` within the removed set at `x`, between the
    // roots there of every boundary the set is read against.
    let column = |w: f64, x: f64| {
        let mut ys: Vec<f64> = vec![0.0, 9.0, 3.0];
        let mut line = |p: (f64, f64), dir: (f64, f64)| {
            if dir.0.abs() > 1e-15 {
                ys.push(p.1 + dir.1 * (x - p.0) / dir.0);
            }
        };
        for off in [0.0, w] {
            line((corner.0 + into.0 * off, corner.1 + into.1 * off), along);
        }
        for (e, t, n) in &ends {
            line(*e, (-t.1, t.0));
            for off in [0.0, w] {
                line((e.0 - n.0 * off, e.1 - n.1 * off), *t);
            }
        }
        for th in [tf, tc] {
            line((cx, cy), (th.cos(), th.sin()));
        }
        if pivot || sector {
            line(corner, into);
            let (c, s) = ((0.5 * turn).cos(), (0.5 * turn).sin());
            for sign in [1.0, -1.0] {
                line(
                    corner,
                    (
                        back.0 * c - sign * back.1 * s,
                        sign * back.0 * s + back.1 * c,
                    ),
                );
            }
            for disc in [w, r] {
                let h = disc * disc - (x - corner.0) * (x - corner.0);
                if h >= 0.0 {
                    ys.push(corner.1 + h.sqrt());
                    ys.push(corner.1 - h.sqrt());
                }
            }
        }
        for k in 0..outline.len() {
            let (a, b) = (outline[k], outline[(k + 1) % outline.len()]);
            if (a.0 - x) * (b.0 - x) <= 0.0 && (b.0 - a.0).abs() > 1e-15 {
                ys.push(a.1 + (b.1 - a.1) * (x - a.0) / (b.0 - a.0));
            }
        }
        for rad in [big, big - w] {
            let h = rad * rad - (x - cx) * (x - cx);
            if h >= 0.0 {
                ys.push(cy + h.sqrt());
                ys.push(cy - h.sqrt());
            }
        }
        ys.retain(|y| (0.0..=9.0).contains(y));
        ys.sort_by(f64::total_cmp);
        ys.windows(2)
            .filter(|p| removed_at(w, x, 0.5 * (p[0] + p[1])))
            .map(|p| p[1] - p[0])
            .sum::<f64>()
    };
    let slice = |phi: f64| {
        let w = r * (1.0 - phi.cos());
        swept::adaptive(&|x| column(w, x), 0.0, 8.0) * r * phi.cos()
    };
    swept::adaptive(&slice, 0.0, core::f64::consts::FRAC_PI_2)
}

/// What rounding takes from [`curved_turn`]'s prism within the sector of
/// radius `r` behind the notch's corner (as [`turn_removed_within`] clips
/// it), through the top `r` of the prism, estimated by sampling: a fixed
/// stream of points uniform over that sector, each counted where `prism`
/// holds it and `rounded` does not, each read against its own boundary
/// meshed to a chord of 1e-4. With one standard error of the estimate.
fn sector_sampled(
    model: &ogeom_topo::Model,
    prism: &ogeom_topo::Shape,
    rounded: &ogeom_topo::Shape,
    turn: f64,
    r: f64,
) -> (f64, f64) {
    use ogeom_algo::{Containment, SolidMesh};
    const SAMPLES: u32 = 20_000;
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let (before, after) = (
        SolidMesh::of(model, prism, fine, T).unwrap(),
        SolidMesh::of(model, rounded, fine, T).unwrap(),
    );
    // The way back from the notch: between the arc's end tangent reversed
    // and the straight edge's direction, as the oracle reads it.
    let (cx, cy) = (6.0_f64, 3.4_f64 - 5.2);
    let tc = (3.0 - cy).atan2(4.0 - cx);
    let ahead = (-tc.sin(), tc.cos());
    let (s, c) = turn.sin_cos();
    let along = (
        ahead.0.mul_add(c, ahead.1 * s),
        (-ahead.0).mul_add(s, ahead.1 * c),
    );
    let back = {
        let b = (ahead.0 - along.0, ahead.1 - along.1);
        let l = b.0.hypot(b.1);
        (b.0 / l, b.1 / l)
    };
    let mid = back.1.atan2(back.0);
    // SplitMix64, seeded: the same points every run.
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut unit = || {
        state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^= z >> 31;
        f64::from(u32::try_from(z >> 32).unwrap()) / 4_294_967_296.0
    };
    let mut taken = 0_u32;
    for _ in 0..SAMPLES {
        let reach = r * unit().sqrt();
        let theta = mid + turn * (unit() - 0.5);
        let z = 2.0 - r * unit();
        let p = Point::new(4.0 + reach * theta.cos(), 3.0 + reach * theta.sin(), z);
        let held = |mesh: &SolidMesh| mesh.holds(p, T).unwrap() == Containment::In;
        if held(&before) && !held(&after) {
            taken += 1;
        }
    }
    let region = 0.5 * turn * r * r * r;
    let share = f64::from(taken) / f64::from(SAMPLES);
    (
        region * share,
        region * (share * (1.0 - share) / f64::from(SAMPLES)).sqrt(),
    )
}

/// The near-closed notch's bounds stand far apart, the sector its rolled
/// corner stands in being much larger than what the corner takes there
/// ([`rounds_at_a_curved_crease_turning_nearly_back_on_itself`]), so what
/// rounding takes in that sector is also sampled, the points read off the
/// prism's and the rounded solid's own boundaries, and held to the exact
/// figure clipped to the sector within four standard errors.
#[test]
#[ignore = "heavy"]
fn a_near_closed_notch_s_rolled_corner_takes_what_sampling_finds() {
    for degrees in [170.0_f64, 178.0] {
        let mut model = ogeom_topo::Model::new();
        let (prism, edges) = curved_turn(&mut model, degrees.to_radians());
        let rounded = ogeom_fillet::fillet_edges(&mut model, &prism, &edges, 0.5, T).unwrap();
        let within = turn_removed_within(degrees.to_radians(), 0.5, Corner::Rolled, true);
        let sampled = sector_sampled(&model, &prism, &rounded.shape, degrees.to_radians(), 0.5);
        println!(
            "{degrees} degrees: in the corner's sector, exact {within:.6}, sampled {:.6} \
             +- {:.6} (one standard error)",
            sampled.0, sampled.1
        );
        assert!(
            (sampled.0 - within).abs() < 4.0 * sampled.1,
            "{degrees} degrees: sampled {sampled:?} against {within}"
        );
    }
}

/// Two top creases, an arc's and a straight edge's, rounded where the
/// outline turns back on itself past them, one request in either order.
/// At 150 degrees the two contacts on the top, carried on past the corner,
/// cross within the solid: the two bands, carried on through the inside
/// corner, take what that rule takes ([`turn_removed_by`] with
/// [`Corner::Joined`]). Turned 170 or 178
/// degrees back they cross only past the solid's far side, and the ball
/// rolls round the notch's edge between the two bands' ends, touching it
/// and the top ([`turn_removed_by`]; ADR 0001, amendment "corners where a
/// round's rail misses"): carried on, the bands' overlap past the corner
/// would meet the notch on the corner's line and pinch the solid. One
/// request is one rounding, so the order the edges are asked in changes
/// nothing, to the bit. Each figure holds the arc's band running on past
/// its far end at (8, 3), where the arc meets the wall x = 8 obliquely,
/// through the corner its ball still finds there.
#[test]
fn rounds_at_a_curved_crease_turning_nearly_back_on_itself() {
    for (degrees, corner) in [
        (150.0_f64, Corner::Joined),
        (170.0, Corner::Rolled),
        (178.0, Corner::Rolled),
    ] {
        let mut model = ogeom_topo::Model::new();
        let (prism, edges) = curved_turn(&mut model, degrees.to_radians());
        let removed = removed_both_ways(&model, &prism, &edges, 0.5);
        let want = turn_removed_by(degrees.to_radians(), 0.5, corner);
        assert_eq!(
            removed[0].to_bits(),
            removed[1].to_bits(),
            "{degrees} degrees: {removed:?}"
        );
        assert!(
            (removed[0] - want).abs() < 2e-6,
            "{degrees} degrees: {removed:?} against {want}"
        );
        if corner != Corner::Rolled {
            continue;
        }
        // Between nothing taken behind the corner and the whole sector to
        // the ball's touch on the edge.
        let (flat, hollow) = (
            turn_removed_by(degrees.to_radians(), 0.5, Corner::Flat),
            turn_removed_by(degrees.to_radians(), 0.5, Corner::Hollow),
        );
        println!(
            "{degrees} degrees: flat {flat:.6}, rolled {:.6} (oracle {want:.6}), hollow \
             {hollow:.6}, gap {:.6}",
            removed[0],
            hollow - flat
        );
        assert!(
            flat < removed[0] && removed[0] < hollow,
            "{degrees} degrees: {flat} < {} < {hollow}",
            removed[0]
        );
        let rounded = ogeom_fillet::fillet_edges(&mut model, &prism, &edges, 0.5, T).unwrap();
        // The rolled corner, the faces both rounds generate, tangent to each
        // band and to the top where it meets them.
        let (a, b) = (
            rounded.history.generated(&edges[0]),
            rounded.history.generated(&edges[1]),
        );
        let rolled: Vec<_> = a
            .iter()
            .filter(|f| b.iter().any(|g| g.is_same(f)))
            .collect();
        assert!(!rolled.is_empty(), "{degrees} degrees: no rolled corner");
        for face in rolled {
            for c in ogeom_fillet::analyse_blend(&model, &rounded.shape, face, 17, T).unwrap() {
                assert!(
                    c.tangency_error < 1e-6 && c.gap < 1e-6,
                    "{degrees} degrees: {c:?}"
                );
            }
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

/// What an inside corner of [`turned_corner`], its second crease leaving at
/// `plan` from the first and its walls leaning by `leans`, takes beyond its
/// two creases' sections along their lengths: within a box round the
/// corner, a column at a time, what lies above both bands carried on less
/// each crease's section along its own edge, past the vertex along neither.
fn leaning_corner_gain(plan: f64, leans: [f64; 2], r: f64) -> f64 {
    let half = 1.6;
    let (s, c) = plan.sin_cos();
    let first = Section::new(0.0, (leans[0].cos(), -leans[0].sin()), r);
    let second = Section::new(0.0, (leans[1].cos(), -leans[1].sin()), r);
    integral((-half, half), (-half, half), 2e3, |x, y| {
        let (ua, ub) = (y, x * s - y * c);
        let [a0, a1] = first.column(ua);
        let [b0, b1] = second.column(ub);
        let (pa, pb) = (first.past(ua, 0.0), second.past(ub, 0.0));
        let above = meet(&[a0, a1, pa], &[b0, b1, pb]);
        let runs: Vec<(f64, f64)> = above
            .iter()
            .map(|&(p, q)| (p.max(0.0), q.min(2.0)))
            .collect();
        let both = span(&runs) - overlap(&runs, &meet(&[pa], &[pb]));
        let own_a = if x >= 0.0 { span(&[a0, a1]) } else { 0.0 };
        let own_b = if x.mul_add(c, y * s) >= 0.0 {
            span(&[b0, b1])
        } else {
            0.0
        };
        both - own_a - own_b
    })
}

/// A block eight square and two high with a pocket one deep inside the
/// convex polygon `outline`, taken anticlockwise, each wall leaning `d`
/// per unit of depth in turn from the side leaving `outline[0]` (a
/// positive draft narrowing the pocket's floor): the block and the
/// pocket's rim creases in that order.
fn polygon_pocket(
    model: &mut ogeom_topo::Model,
    outline: &[(f64, f64)],
    d: &[f64],
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
    let n = outline.len();
    #[allow(clippy::cast_precision_loss, reason = "a corner count")]
    let share = 1.0 / n as f64;
    let (mut cx, mut cy) = (0.0, 0.0);
    let mut middles = Vec::with_capacity(n);
    let mut sides = Vec::with_capacity(n + 1);
    for k in 0..n {
        let (p, q) = (outline[k], outline[(k + 1) % n]);
        let (dx, dy) = (q.0 - p.0, q.1 - p.1);
        let length = dx.hypot(dy);
        let middle = Point::new(f64::midpoint(p.0, q.0), f64::midpoint(p.1, q.1), 2.0);
        middles.push(middle);
        cx += p.0 * share;
        cy += p.1 * share;
        // The wall's normal into the pocket, tipped by its draft.
        sides.push(past(
            model,
            middle,
            Vector::new(-dy / length, dx / length, d[k]),
        ));
    }
    sides.push(past(model, Point::new(cx, cy, 1.0), Vector::Z));
    let mut pocket = sides[0].clone();
    for side in &sides[1..] {
        pocket = ogeom_bool::common(model, &pocket, side, T).unwrap().shape;
    }
    let cut = ogeom_bool::cut(model, &block, &pocket, T).unwrap().shape;
    let block = explore_unique(model, &cut, ShapeType::Solid).unwrap()[0].clone();
    let edges = middles
        .iter()
        .map(|&at| placed_edge_through(model, &block, at))
        .collect();
    (block, edges)
}

/// What rounding every rim crease of [`polygon_pocket`] at `r` takes: each
/// crease's section along its length, and each corner what its two bands
/// take beyond those sections.
fn pocket_rounding(outline: &[(f64, f64)], d: &[f64], r: f64) -> f64 {
    let lean = |d: f64| -d.atan();
    let n = outline.len();
    let mut want = 0.0;
    for k in 0..n {
        let (p, q) = (outline[k], outline[(k + 1) % n]);
        let area = Section::new(0.0, (lean(d[k]).cos(), -lean(d[k]).sin()), r).area;
        want += (q.0 - p.0).hypot(q.1 - p.1) * area;
        // The corner at `p`: its first crease leaves along the side to `q`,
        // the material on its right, and its second along the side in.
        let before = outline[(k + n - 1) % n];
        let plan = ((before.1 - p.1).atan2(before.0 - p.0) - (q.1 - p.1).atan2(q.0 - p.0))
            .rem_euclid(core::f64::consts::TAU);
        want += leaning_corner_gain(plan, [lean(d[k]), lean(d[(k + n - 1) % n])], r);
    }
    want
}

/// Pockets whose walls lean by different drafts rounded round their
/// rims: four walls each its own way, square in plan, as a trapezoid with
/// corners at 63 and 117 degrees, and as a pentagon; and a pocket a
/// millimetre and a half wide whose two long walls lean. Every corner is
/// crossed or flush, each band crossed at both ends cut in two, and each
/// crease takes its section along its length and each corner what lies
/// above both its bands beyond those sections, in either order. Off square,
/// the ring's wedge holds faces at the acute corners narrower than the gap
/// its halves are sewn across, which the sewing took for slits, and the
/// boolean read a short fitted edge along the pocket's corner as a section
/// of its own: both refused, and every corner was left capped.
#[test]
fn pockets_whose_walls_lean_their_own_ways_meet_at_every_corner() {
    // Each pocket's outline and its walls' drafts.
    type Pocket<'a> = (&'a [(f64, f64)], &'a [f64]);
    let cases: [Pocket; 4] = [
        (
            &[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)],
            &[0.1, 0.05, -0.05, 0.08],
        ),
        (
            &[(1.5, 2.0), (6.5, 2.0), (5.0, 5.0), (3.0, 5.0)],
            &[0.1, 0.05, -0.05, 0.08],
        ),
        (
            &[(2.0, 2.0), (6.0, 2.0), (7.0, 4.5), (4.0, 6.5), (1.0, 4.5)],
            &[0.1, 0.05, -0.05, 0.08, 0.0],
        ),
        (
            &[(2.0, 3.25), (6.0, 3.25), (6.0, 4.75), (2.0, 4.75)],
            &[0.1, 0.0, 0.05, 0.0],
        ),
    ];
    for (outline, d) in cases {
        let mut model = ogeom_topo::Model::new();
        let (block, edges) = polygon_pocket(&mut model, outline, d);
        let want = pocket_rounding(outline, d, 0.5);
        // Each corner's box is integrated to about a micron cubed.
        for got in removed_both_ways(&model, &block, &edges, 0.5) {
            assert!((got - want).abs() < 1e-5, "{d:?}: {got} against {want}");
        }
    }
}

/// Moulded pockets, every wall drafted alike so its foot runs out into the
/// pocket's floor, rounded round their rims: square, and a trapezoid with
/// acute and obtuse corners; and a square one with two adjacent walls
/// drafted alike. Where two walls lean alike, their bands' end sections
/// reach past the plane bisecting the corner, in which they meet flush;
/// each crease takes its section along its length and each corner what
/// lies above both its bands beyond those sections, in either order. The
/// bands could not reach back to that plane from their caps, and every
/// such corner was left capped.
#[test]
fn moulded_pockets_meet_at_every_corner() {
    let square = [(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)];
    let trapezoid = [(1.5, 2.0), (6.5, 2.0), (5.0, 5.0), (3.0, 5.0)];
    for (outline, d) in [
        (square, [0.1; 4]),
        (trapezoid, [0.1; 4]),
        (square, [0.1, 0.1, 0.0, 0.0]),
    ] {
        let mut model = ogeom_topo::Model::new();
        let (block, edges) = polygon_pocket(&mut model, &outline, &d);
        let want = pocket_rounding(&outline, &d, 0.5);
        for got in removed_both_ways(&model, &block, &edges, 0.5) {
            assert!(
                (got - want).abs() < 1e-5,
                "{outline:?} {d:?}: {got} against {want}"
            );
        }
    }
}

/// The top creases of `solid` at `z = 2` leaving `corner`, the circular
/// one first.
fn creases_at(
    model: &ogeom_topo::Model,
    solid: &ogeom_topo::Shape,
    corner: Point,
) -> [ogeom_topo::Shape; 2] {
    let point = |v: &ogeom_topo::Shape| {
        let p = model.node(v).unwrap().data().as_vertex().unwrap().point;
        v.transform(model.datums()).unwrap().apply(p)
    };
    let mut edges: Vec<ogeom_topo::Shape> = explore_unique(model, solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(model, e).unwrap().unwrap();
            let (p, q) = (point(&a), point(&b));
            (p.z - 2.0).abs() < 1e-9
                && (q.z - 2.0).abs() < 1e-9
                && (p.distance(corner) < 1e-7 || q.distance(corner) < 1e-7)
        })
        .collect();
    edges.sort_by_key(|e| {
        let (curve, _) = edge_curve_of(model, e);
        !matches!(curve, ogeom_geom::Curve::Circle(_))
    });
    let [first, second] = edges.as_slice() else {
        panic!("two top creases at the corner");
    };
    [first.clone(), second.clone()]
}

/// The box past the plane through `on` with unit `normal`, forty a side.
fn past_plane(model: &mut ogeom_topo::Model, on: Point, normal: Vector) -> ogeom_topo::Shape {
    let seed = if normal.x.abs() < 0.9 {
        Vector::X
    } else {
        Vector::Y
    };
    let x = seed - normal * seed.dot(normal);
    let x = x / x.magnitude();
    let y = normal.cross(x);
    let frame = Frame::new(
        on - x * 20.0 - y * 20.0,
        Direction::new(normal, T).unwrap(),
        Direction::new(x, T).unwrap(),
        T,
    )
    .unwrap();
    ogeom_algo::make_box(model, frame, (40.0, 40.0, 40.0), T)
        .unwrap()
        .shape
}

/// A drum of radius 3 at its top about the `z` axis, two high, drafted by
/// `k` over the height so its foot stands out for a positive `k` (a cone;
/// upright, a cylinder), fused with a block `x <= 0` whose top crease `y =
/// b` meets the drum's rim at `(-sqrt(9 - b^2), b, 2)`, its wall leaning by
/// `lean` from the vertical, its foot in under the block for a positive
/// lean, the two tops unified into one face. Returns the solid, that
/// corner, and the rim's and the block's top creases leaving it.
fn drum_on_a_leaning_block(
    model: &mut ogeom_topo::Model,
    b: f64,
    lean: f64,
    k: f64,
) -> (ogeom_topo::Shape, Point, [ogeom_topo::Shape; 2]) {
    let frame = Frame::new(Point::ORIGIN, Direction::Z, -Direction::X, T).unwrap();
    let drum = if k == 0.0 {
        ogeom_algo::make_cylinder(model, frame, 3.0, 2.0, T)
    } else {
        ogeom_algo::make_cone(model, frame, 2.0_f64.mul_add(k, 3.0), 3.0, 2.0, T)
    }
    .unwrap()
    .shape;
    let block = ogeom_algo::make_box(
        model,
        Frame::new(Point::new(-7.0, -5.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
        (7.0, 6.0 + b, 2.0),
        T,
    )
    .unwrap()
    .shape;
    let past = past_plane(
        model,
        Point::new(0.0, b, 2.0),
        Vector::new(0.0, lean.cos(), -lean.sin()),
    );
    let block = ogeom_bool::cut(model, &block, &past, T).unwrap().shape;
    let fused = ogeom_bool::fuse(model, &drum, &block, T).unwrap().shape;
    let fused = explore_unique(model, &fused, ShapeType::Solid).unwrap()[0].clone();
    let unified = ogeom_heal::unify_same_domain(model, &fused, T)
        .unwrap()
        .0
        .shape;
    let solid = explore_unique(model, &unified, ShapeType::Solid).unwrap()[0].clone();
    let corner = Point::new(-(9.0 - b * b).sqrt(), b, 2.0);
    let edges = creases_at(model, &solid, corner);
    (solid, corner, edges)
}

/// A block eight square and two high with a pocket one deep: the disc of
/// radius 2.5 about `(4, 4)` up to its chord `y = 4 + c` at the top, that
/// wall drafted by `d` over the depth so the floor narrows. Returns the
/// block, the rim's corner `(4 + sqrt(2.5^2 - c^2), 4 + c, 2)`, and the
/// arc's and the chord's rim creases leaving it.
fn drafted_d_pocket(
    model: &mut ogeom_topo::Model,
    c: f64,
    d: f64,
) -> (ogeom_topo::Shape, Point, [ogeom_topo::Shape; 2]) {
    let big = 2.5;
    let block = ogeom_algo::make_box(model, Frame::WORLD, (8.0, 8.0, 2.0), T)
        .unwrap()
        .shape;
    let disc = ogeom_algo::make_cylinder(
        model,
        Frame::new(Point::new(4.0, 4.0, 1.0), Direction::Z, Direction::Y, T).unwrap(),
        big,
        2.0,
        T,
    )
    .unwrap()
    .shape;
    let side = past_plane(
        model,
        Point::new(4.0, 4.0 + c, 2.0),
        Vector::new(0.0, -1.0, d),
    );
    let pocket = ogeom_bool::common(model, &disc, &side, T).unwrap().shape;
    let cut = ogeom_bool::cut(model, &block, &pocket, T).unwrap().shape;
    let solid = explore_unique(model, &cut, ShapeType::Solid).unwrap()[0].clone();
    let corner = Point::new(4.0 + big.mul_add(big, -c * c).sqrt(), 4.0 + c, 2.0);
    let edges = creases_at(model, &solid, corner);
    (solid, corner, edges)
}

/// A corner where a top crease on a drum of radius `big` about `centre`
/// meets a straight one at `at`, both under the top `z = 2`: the material
/// within the drum where `within`, the drum's wall drafted by `k` (its
/// outward normal in the half-plane through the axis `(1, k)`, off the
/// axis and up; upright for zero), the drum's crease carried on past the
/// vertex along `ahead`, the straight wall's outward normal in plan `side`
/// and its lean `lean`. Within the square of half side `half`
/// round the corner, what lies above both bands carried on and in the
/// material, a column at a time: the drum's band its torus band before its
/// cap and the straight band of its end section past it, out of the open
/// corner off the drum and past the straight wall.
#[allow(clippy::too_many_arguments, reason = "one corner, all its data")]
fn turned_corner_removed(
    at: (f64, f64),
    (centre, big, within): ((f64, f64), f64, bool),
    ahead: (f64, f64),
    side: (f64, f64),
    (lean, k): (f64, f64),
    r: f64,
    half: f64,
) -> f64 {
    let towards = if within { 1.0 } else { -1.0 };
    let normal = (
        towards * (at.0 - centre.0) / big,
        towards * (at.1 - centre.1) / big,
    );
    let turned = Section::new(0.0, (1.0, k), r);
    let straight = Section::new(0.0, (lean.cos(), -lean.sin()), r);
    integral(
        (at.0 - half, at.0 + half),
        (at.1 - half, at.1 + half),
        2e3,
        |x, y| {
            let (px, py) = (x - at.0, y - at.1);
            let off = towards * ((x - centre.0).hypot(y - centre.1) - big);
            let ua = if px.mul_add(ahead.0, py * ahead.1) > 0.0 {
                px.mul_add(normal.0, py * normal.1)
            } else {
                off
            };
            let ub = px.mul_add(side.0, py * side.1);
            let [a0, a1] = turned.column(ua);
            let [b0, b1] = straight.column(ub);
            let (pa, pb) = (turned.past(ua, 0.0), straight.past(ub, 0.0));
            let above = meet(&[a0, a1, pa], &[b0, b1, pb]);
            let runs: Vec<(f64, f64)> = above
                .iter()
                .map(|&(p, q)| (p.max(0.0), q.min(2.0)))
                .collect();
            let open = turned.past(off, 0.0);
            span(&runs) - overlap(&runs, &meet(&[open], &[pb]))
        },
    )
}

/// What rounding `edges` of `solid` at `r` takes within the square of half
/// side `half` round `at`, in the order given and reversed, each result
/// checked valid.
fn removed_near(
    model: &ogeom_topo::Model,
    solid: &ogeom_topo::Shape,
    edges: &[ogeom_topo::Shape; 2],
    r: f64,
    at: Point,
    half: f64,
) -> [f64; 2] {
    [edges.clone(), [edges[1].clone(), edges[0].clone()]].map(|named| {
        let mut model = model.clone();
        let rounded = ogeom_fillet::fillet_edges(&mut model, solid, &named, r, T).unwrap();
        assert!(
            ogeom_algo::check(&model, &rounded.shape, T)
                .unwrap()
                .is_valid()
        );
        let window = ogeom_algo::make_box(
            &mut model,
            Frame::new(
                Point::new(at.x - half, at.y - half, -1.0),
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
        let before = ogeom_bool::common(&mut model, solid, &window, T)
            .unwrap()
            .shape;
        let after = ogeom_bool::common(&mut model, &rounded.shape, &window, T)
            .unwrap()
            .shape;
        volume(&model, &before, 1e-4) - volume(&model, &after, 1e-4)
    })
}

/// Two top creases rounded where a drum fused onto a block meets the
/// block's leaning wall: the rim's crease turns about the drum's axis, the
/// block's is straight, and their end sections differ, so the junction is
/// crossed. The rim's band is the torus band of its section before its cap
/// and, carried on, the straight band of its end section past it; the
/// corner takes what lies above both bands in the material, measured
/// within a box round the corner a column at a time, in either order.
/// Without the junction the two bands met as each was built, and left the
/// corner standing.
#[test]
fn rounds_at_a_drum_meeting_a_leaning_wall_take_all_above_both() {
    let r = 0.5;
    let half = 1.2;
    for (b, lean) in [(1.0_f64, 0.1_f64.atan()), (2.0, 0.1_f64.atan())] {
        let mut model = ogeom_topo::Model::new();
        let (solid, at, edges) = drum_on_a_leaning_block(&mut model, b, lean, 0.0);
        let reach = (9.0 - b * b).sqrt();
        let want = turned_corner_removed(
            (at.x, at.y),
            ((0.0, 0.0), 3.0, true),
            (-b / 3.0, -reach / 3.0),
            (0.0, 1.0),
            (lean, 0.0),
            r,
            half,
        );
        for got in removed_near(&model, &solid, &edges, r, at, half) {
            assert!(
                (got - want).abs() < 5e-6,
                "b {b}, leaning {lean}: {got} against {want}"
            );
        }
    }
}

/// As [`rounds_at_a_drum_meeting_a_leaning_wall_take_all_above_both`], the
/// block's wall leaning out from under it. The straight band carried on
/// then stands off the drum before the rim's cap, under the block's top,
/// where the material lies above both bands: the rim's band widened past
/// the drum takes it, a widening on the drum's own cylinder. Without it
/// the request fell back to capped bands and left that material standing.
/// Its junction's booleans build without the boolean's last retry, which
/// ends contacts where their curves part by their own radii.
#[test]
fn rounds_at_a_drum_meeting_a_wall_leaning_out_take_all_above_both() {
    let r = 0.5;
    let half = 1.2;
    let lean = -0.1_f64.atan();
    for b in [1.0_f64, 2.0] {
        let mut model = ogeom_topo::Model::new();
        let (solid, at, edges) = drum_on_a_leaning_block(&mut model, b, lean, 0.0);
        let reach = (9.0 - b * b).sqrt();
        let want = turned_corner_removed(
            (at.x, at.y),
            ((0.0, 0.0), 3.0, true),
            (-b / 3.0, -reach / 3.0),
            (0.0, 1.0),
            (lean, 0.0),
            r,
            half,
        );
        let (removed, work) =
            ogeom_bool::work::counted(|| removed_near(&model, &solid, &edges, r, at, half));
        assert_eq!(work.stated_retries, 0, "b {b}: the last retry was entered");
        for got in removed {
            assert!(
                (got - want).abs() < 5e-6,
                "b {b}, leaning {lean}: {got} against {want}"
            );
        }
    }
}

/// As [`rounds_at_a_drum_meeting_a_leaning_wall_take_all_above_both`], the
/// drum drafted so it narrows to its foot (a cone), the block's wall
/// square or leaning either way. The rim's round is a torus band, its
/// section the same at every turn about the axis, so the junction takes
/// what lies above both bands, the open corner bounded by the cone itself.
/// The rim's round was fitted, and the junction took no fitted band: the
/// request refused in both orders, or fell 3.6e-3 to 1.31e-2 short.
#[test]
fn rounds_at_a_drafted_drum_meeting_a_wall_take_all_above_both() {
    let r = 0.5;
    let half = 1.2;
    let k = -0.1;
    for b in [1.0_f64, 2.0] {
        for lean in [0.0, 0.1_f64.atan(), -0.1_f64.atan()] {
            let mut model = ogeom_topo::Model::new();
            let (solid, at, edges) = drum_on_a_leaning_block(&mut model, b, lean, k);
            let reach = (9.0 - b * b).sqrt();
            let want = turned_corner_removed(
                (at.x, at.y),
                ((0.0, 0.0), 3.0, true),
                (-b / 3.0, -reach / 3.0),
                (0.0, 1.0),
                (lean, k),
                r,
                half,
            );
            for got in removed_near(&model, &solid, &edges, r, at, half) {
                assert!(
                    (got - want).abs() < 5e-6,
                    "b {b}, leaning {lean}: {got} against {want}"
                );
            }
        }
    }
}

/// As [`rounds_at_a_drafted_drum_meeting_a_wall_take_all_above_both`], the
/// drum drafted so its foot stands out. With the block's wall leaning in,
/// and at the farther crease square, the junction's reach was read off the
/// whole conic that the cone meets the block's wall plane in rather than
/// off the stretch the corner's edge uses: the straight band was carried on
/// across the drum and out of its far side, and the corner refused. With
/// the wall leaning out too, the two walls lean alike with their feet out,
/// the rim's end section stands past the plane bisecting the turn, and the
/// flush meeting in that plane refused. Each fell back to the bands as
/// built, 3.5e-3 to 7.5e-3 short. At the nearer crease with the wall
/// square the band carried on past its cap crosses the cone at a shallow
/// turn beside the wall. Swept from the cap's arc spelled as a spline,
/// which may be a fit, that crossing was read as a touch; swept from the
/// cap's own circle, the march found no seed in the window it walks the
/// band in. The corner refused, and fell 6.5e-3 short.
#[test]
fn rounds_at_a_drum_drafted_foot_out_meeting_a_wall_take_all_above_both() {
    let r = 0.5;
    let half = 1.2;
    let k = 0.1;
    let (inward, outward) = (0.1_f64.atan(), -0.1_f64.atan());
    for (b, lean) in [
        (1.0_f64, 0.0),
        (1.0, inward),
        (1.0, outward),
        (2.0, 0.0),
        (2.0, inward),
        (2.0, outward),
    ] {
        let mut model = ogeom_topo::Model::new();
        let (solid, at, edges) = drum_on_a_leaning_block(&mut model, b, lean, k);
        let reach = (9.0 - b * b).sqrt();
        let want = turned_corner_removed(
            (at.x, at.y),
            ((0.0, 0.0), 3.0, true),
            (-b / 3.0, -reach / 3.0),
            (0.0, 1.0),
            (lean, k),
            r,
            half,
        );
        for got in removed_near(&model, &solid, &edges, r, at, half) {
            assert!(
                (got - want).abs() < 5e-6,
                "b {b}, leaning {lean}: {got} against {want}"
            );
        }
    }
}

/// A drum drafted with its foot out (a cone) fused onto a block whose wall
/// leans out, its two top creases rounded one call each, in either order.
/// Separate rounds depend on their order, as in every application (ADR
/// 0001, amendment "corners where a round's rail misses"), so each order is
/// read against its own rule, a column at a time within a box round the
/// corner ([`cone_rim_corner_removed`]):
/// - the block's crease alone ends flush on the drum's cone (amendment
///   "M5e", a free end into the material), its section kept outside it;
/// - the rim alone ends flush on the block's leaning wall there, its band
///   kept on the wall's open side near the corner;
/// - the rim after the crease runs on past the corner round its circle
///   while the crease's cut leaves its contact on the top in the solid and
///   its contact on the cone carried on outside it, and ends in the plane
///   through the axis there (amendment "corners where a round's rail
///   misses");
/// - the crease after the rim runs on along the top and the piece of the
///   wall the rim's flush end left, one edge, until the rim's band trims
///   it (amendment "M5e"; #1098's rule for a band ending at an earlier
///   one), its section kept short of the band's tube.
///
/// The two orders differ only within the box. The current path took one
/// corner in either order, its own figure. The crease after the rim was
/// built by the current path, its own figure, while the rim's flush end
/// left the crease's line split at the old corner.
#[test]
fn a_cone_rim_that_meets_unresolved_is_rounded_fitted() {
    let r = 0.5;
    let half = 1.2;
    let mut model = ogeom_topo::Model::new();
    let (solid, at, edges) = drum_on_a_leaning_block(&mut model, 2.0, -0.1_f64.atan(), 0.1);
    let window = |model: &mut ogeom_topo::Model| {
        ogeom_algo::make_box(
            model,
            Frame::new(
                Point::new(at.x - half, at.y - half, -1.0),
                Direction::Z,
                Direction::X,
                T,
            )
            .unwrap(),
            (half * 2.0, half * 2.0, 4.0),
            T,
        )
        .unwrap()
        .shape
    };
    let within = |model: &mut ogeom_topo::Model, shape: &ogeom_topo::Shape| {
        let window = window(model);
        let kept = ogeom_bool::common(model, shape, &window, T).unwrap().shape;
        volume(model, &kept, 1e-4)
    };
    let mut scratch = model.clone();
    let before = within(&mut scratch, &solid);
    let total_before = volume(&model, &solid, 1e-4);
    let want = |rule: ConeRimCorner| cone_rim_corner_removed((at.x, at.y), half, r, rule);
    let close = |what: &str, got: f64, want: f64| {
        assert!((got - want).abs() < 1e-6, "{what}: {got} against {want}");
    };
    let mut totals = Vec::new();
    let mut windows = Vec::new();
    for order in [[1, 0], [0, 1]] {
        let mut model = model.clone();
        let once = ogeom_fillet::fillet_edge(&mut model, &solid, &edges[order[0]], r, T).unwrap();
        let first = before - within(&mut model, &once.shape);
        if order[0] == 1 {
            close("the crease alone", first, want(ConeRimCorner::CreaseAlone));
        } else {
            close("the rim alone", first, want(ConeRimCorner::RimAlone));
        }
        let [later] = once.history.trace(&edges[order[1]]) else {
            panic!("the later crease is not one edge after the first round");
        };
        let later = later.clone();
        let rounded = ogeom_fillet::fillet_edge(&mut model, &once.shape, &later, r, T).unwrap();
        assert!(
            ogeom_algo::check(&model, &rounded.shape, T)
                .unwrap()
                .is_valid()
        );
        let both = before - within(&mut model, &rounded.shape);
        if order[0] == 1 {
            close(
                "the rim after the crease",
                both,
                want(ConeRimCorner::RimAfterCrease),
            );
        } else {
            close(
                "the crease after the rim",
                both,
                want(ConeRimCorner::CreaseAfterRim),
            );
        }
        windows.push(both);
        totals.push(total_before - volume(&model, &rounded.shape, 1e-4));
    }
    // Past the box the two orders take the same.
    let outside = (totals[0] - windows[0]) - (totals[1] - windows[1]);
    assert!(
        outside.abs() < 1e-6,
        "the orders differ past the corner by {outside}"
    );
}

/// Which step of [`a_cone_rim_that_meets_unresolved_is_rounded_fitted`]
/// [`cone_rim_corner_removed`] reads.
#[derive(Clone, Copy)]
enum ConeRimCorner {
    CreaseAlone,
    RimAlone,
    RimAfterCrease,
    CreaseAfterRim,
}

/// What a step of [`a_cone_rim_that_meets_unresolved_is_rounded_fitted`]
/// takes within the square of half side `half` round the corner `at`, a
/// column at a time. The drum's radius is `3 + 0.1 (2 - z)`; the block's
/// wall through the crease `y = 2` leans by `atan 0.1` with its foot out.
/// The crease's section is the round's in the plane square to it, kept
/// outside the drum (where its column stands above `2 - 10 (rho - 3)`);
/// the rim's is its torus band's, over the arc from the corner round the
/// drum (polar angle at most the corner's), kept on the wall's open side
/// when the rim is rounded alone, and run on past the corner after the
/// crease to where it stops. After the rim, the crease runs on along the
/// top and the wall's piece the rim left, to the rim's torus band, and
/// ends flush on it: a point of its section is taken while the line
/// through it along the crease has not yet reached the band's tube (the
/// tube's centre circle of radius `3 + u`, `u` the rim section's centre,
/// at the height of that centre, its radius `r`).
fn cone_rim_corner_removed(at: (f64, f64), half: f64, r: f64, rule: ConeRimCorner) -> f64 {
    let k = 0.1;
    let lean = -0.1_f64.atan();
    let (ny, nz) = (lean.cos(), -lean.sin());
    let crease = Section::new(2.0, (ny, nz), r);
    let rim = Section::new(0.0, (1.0, k), r);
    let corner = at.1.atan2(at.0);
    // After the crease, the rim's band runs on past the corner round its
    // circle while its contact on the top stands in the solid and its
    // contact on the cone carried on stands outside it (amendment "corners
    // where a round's rail misses"), to the angle `stop`.
    let stop = match rule {
        ConeRimCorner::RimAfterCrease => {
            let eps = 1e-9;
            let in_crease = |x: f64, y: f64, z: f64| {
                let rho = x.hypot(y);
                rho > 3.0 + eps
                    && z > 2.0 - (rho - 3.0) / k + eps
                    && crease
                        .column(y)
                        .iter()
                        .any(|&(lo, hi)| z > lo + eps && z < hi - eps)
            };
            let in_solid = |x: f64, y: f64, z: f64| {
                let in_height = (-eps..=2.0 + eps).contains(&z);
                let drum = x.hypot(y) <= 3.0 + k * (2.0 - z) + eps;
                let block = x <= eps && ny.mul_add(y - 2.0, nz * (z - 2.0)) <= eps;
                in_height && (drum || block) && !in_crease(x, y, z)
            };
            let m: (f64, f64) = (1.0, k);
            let m = (m.0 / m.0.hypot(m.1), m.1 / m.0.hypot(m.1));
            let cap = (rim.centre.0, 2.0);
            let wall = (rim.centre.0 + m.0 * r, rim.centre.1 + m.1 * r);
            let trims = |theta: f64| {
                let at = |(u, z): (f64, f64)| ((3.0 + u) * theta.cos(), (3.0 + u) * theta.sin(), z);
                let (cx, cy, cz) = at(cap);
                let (wx, wy, wz) = at(wall);
                in_solid(cx, cy, cz) && !in_solid(wx, wy, wz)
            };
            let step = 0.25 * r / 3.0;
            let (mut good, mut bad) = (corner, corner + step);
            while trims(bad) {
                good = bad;
                bad += step;
            }
            for _ in 0..60 {
                let mid = f64::midpoint(good, bad);
                if trims(mid) {
                    good = mid;
                } else {
                    bad = mid;
                }
            }
            good
        }
        _ => corner,
    };
    integral(
        (at.0 - half, at.0 + half),
        (at.1 - half, at.1 + half),
        2e3,
        |x, y| {
            let rho = x.hypot(y);
            let outside_drum = if rho > 3.0 {
                vec![((2.0 - (rho - 3.0) / k).max(0.0), 2.0)]
            } else {
                Vec::new()
            };
            let crease_runs = meet(&crease.column(y), &outside_drum);
            let on_arc = y.atan2(x) <= stop;
            let rim_runs: Vec<(f64, f64)> = if on_arc {
                rim.column(rho - 3.0).to_vec()
            } else {
                Vec::new()
            };
            match rule {
                ConeRimCorner::CreaseAlone => span(&crease_runs),
                ConeRimCorner::RimAlone => {
                    // Behind the wall's plane: z at most 2 - (ny / nz)(y - 2).
                    let open = [(2.0 - ny / nz * (y - 2.0), f64::INFINITY)];
                    span(&meet(&rim_runs, &open))
                }
                ConeRimCorner::RimAfterCrease => {
                    span(&crease_runs) + span(&rim_runs) - overlap(&crease_runs, &rim_runs)
                }
                ConeRimCorner::CreaseAfterRim => {
                    let open = [(2.0 - ny / nz * (y - 2.0), f64::INFINITY)];
                    let rim_part = meet(&rim_runs, &open);
                    // Before the tube along the crease's line: farther out
                    // than the tube at this height, or at a height the
                    // tube does not reach.
                    let (rc, zc) = (3.0 + rim.centre.0, rim.centre.1);
                    let below = (f64::NEG_INFINITY, zc - r);
                    let above = (zc + r, f64::INFINITY);
                    let before: Vec<(f64, f64)> = if rho - rc >= r {
                        vec![(f64::NEG_INFINITY, f64::INFINITY)]
                    } else if rho > rc {
                        let h = r.mul_add(r, -(rho - rc).powi(2)).sqrt();
                        vec![(f64::NEG_INFINITY, zc - h), (zc + h, f64::INFINITY)]
                    } else {
                        vec![below, above]
                    };
                    let crease_part = meet(&crease.column(y), &before);
                    span(&crease_part) + span(&rim_part) - overlap(&crease_part, &rim_part)
                }
            }
        },
    )
}

/// A pocket of a disc up to a chord whose wall is drafted, the floor
/// narrower: the arc's band turns about the bore's axis, the material
/// outside it, and the chord's band is straight with another section, so
/// the junction is crossed. The open corner past the chord's wall is what
/// stands within the bore, which bends away from its tangent plane into
/// the material, so the corner is cut by the bore's own cylinder. Both
/// orders take what lies above both bands in the material, with the chord
/// off the centre either way and through it. A bore took no junction, and
/// the request fell back to capped bands. With the chord through the centre
/// or past it, the drafted wall leans under the host off the bore's
/// cylinder, within it, beside the chord's own band before the arc's cap,
/// where the arc's band was widened only round its own part: the junction
/// was refused, and the capped bands fell short or refused in the boolean
/// in one naming order or both.
#[test]
fn rounds_at_a_bore_meeting_a_drafted_chord_take_all_above_both() {
    let r = 0.5;
    let half = 1.2;
    for (c, d) in [(0.5_f64, 0.1_f64), (0.0, 0.1), (0.0, 0.05), (-0.5, 0.1)] {
        let mut model = ogeom_topo::Model::new();
        let (solid, at, edges) = drafted_d_pocket(&mut model, c, d);
        let reach = 2.5_f64.mul_add(2.5, -c * c).sqrt();
        let want = turned_corner_removed(
            (at.x, at.y),
            ((4.0, 4.0), 2.5, false),
            (-c / 2.5, reach / 2.5),
            (0.0, -1.0),
            (-d.atan(), 0.0),
            r,
            half,
        );
        for got in removed_near(&model, &solid, &edges, r, at, half) {
            assert!(
                (got - want).abs() < 5e-6,
                "chord at {c}, draft {d}: {got} against {want}"
            );
        }
    }
}

/// A pocket of a disc up to a chord, its two rim creases rounded: the arc
/// and the chord meet at both their ends, a loop of two, each corner flush
/// in its bisecting plane. With the chord off the centre and on it, where
/// the arc's two caps lie in one plane through the axis, each corner takes
/// what lies above both bands carried on, measured within a box round one
/// corner a column at a time, in either order. The loop's second junction
/// was dropped, and its corner kept its caps; on the centre, the cap at the
/// far corner was read for the near one and neither corner met.
#[test]
fn a_pocket_of_a_disc_and_its_chord_meets_at_both_corners() {
    let r = 0.5;
    let half = 1.2;
    for c in [0.5_f64, 0.0] {
        let mut model = ogeom_topo::Model::new();
        let (solid, at, edges) = drafted_d_pocket(&mut model, c, 0.0);
        let reach = 2.5_f64.mul_add(2.5, -c * c).sqrt();
        let want = turned_corner_removed(
            (at.x, at.y),
            ((4.0, 4.0), 2.5, false),
            (-c / 2.5, reach / 2.5),
            (0.0, -1.0),
            (0.0, 0.0),
            r,
            half,
        );
        for got in removed_near(&model, &solid, &edges, r, at, half) {
            assert!(
                (got - want).abs() < 5e-6,
                "chord at {c}: {got} against {want}"
            );
        }
    }
}

/// A disc of radius 2.5 about `(4, 4)` at the top, two high, its drum
/// drafted by `k` over the height so its foot stands out for a positive
/// `k` (a cone; upright, a cylinder), up to its chord `y = 4 + c` at the
/// top, the chord's wall leaning `d` per unit of depth, its foot in for a
/// positive `d`. Returns the block, the rim's corner `(4 + sqrt(2.5^2 -
/// c^2), 4 + c, 2)`, and the arc's and the chord's top creases leaving it.
fn d_shaped_block(
    model: &mut ogeom_topo::Model,
    c: f64,
    d: f64,
    k: f64,
) -> (ogeom_topo::Shape, Point, [ogeom_topo::Shape; 2]) {
    let big = 2.5;
    let frame = Frame::new(Point::new(4.0, 4.0, 0.0), Direction::Z, Direction::Y, T).unwrap();
    let disc = if k == 0.0 {
        ogeom_algo::make_cylinder(model, frame, big, 2.0, T)
    } else {
        ogeom_algo::make_cone(model, frame, 2.0_f64.mul_add(k, big), big, 2.0, T)
    }
    .unwrap()
    .shape;
    let normal = Vector::new(0.0, 1.0, d) * (1.0 / d.hypot(1.0));
    let past = past_plane(model, Point::new(4.0, 4.0 + c, 2.0), normal);
    let cut = ogeom_bool::cut(model, &disc, &past, T).unwrap().shape;
    let solid = explore_unique(model, &cut, ShapeType::Solid).unwrap()[0].clone();
    let corner = Point::new(4.0 + big.mul_add(big, -c * c).sqrt(), 4.0 + c, 2.0);
    let edges = creases_at(model, &solid, corner);
    (solid, corner, edges)
}

/// Where the arc's band of [`d_shaped_block`] stops past the vertex, alone
/// (ADR 0001, amendment "corners where a round's rail misses"): it runs on
/// round the axis while its contact on the top stands in the solid, and
/// ends in the plane through the axis where that contact, a circle of
/// radius `big + top` about `(4, 4)` (`top` the contact's offset from the
/// rim, negative), crosses the chord's top line `y = 4 + c`. The plane's
/// unit normal pointing on past it: a point `(x, y)` is past the stop
/// where `((x - 4), (y - 4))` dotted with it is positive.
fn arc_stop(c: f64, big: f64, top: f64) -> (f64, f64) {
    let rho = big + top;
    let x = rho.mul_add(rho, -c * c).sqrt();
    (-c / rho, x / rho)
}

/// What rounding the arc and the chord of [`d_shaped_block`] takes within
/// the square of half side `half` round its corner `at`, a column at a
/// time: the union of what each band takes as it does alone, in the
/// material. The arc's band is its torus band run on past the vertex to
/// [`arc_stop`]; the chord's is its straight band up to
/// its own cap, square to the chord at the vertex.
fn d_shaped_corner_removed(at: (f64, f64), c: f64, r: f64, half: f64) -> f64 {
    let (centre, big) = ((4.0, 4.0), 2.5);
    let section = Section::new(0.0, (1.0, 0.0), r);
    let ahead = arc_stop(c, big, section.centre.0);
    integral(
        (at.0 - half, at.0 + half),
        (at.1 - half, at.1 + half),
        2e3,
        |x, y| {
            let off = (x - centre.0).hypot(y - centre.1) - big;
            if off > 0.0 || y > at.1 {
                return 0.0;
            }
            let (px, py) = (x - at.0, y - at.1);
            let arc = if (x - centre.0).mul_add(ahead.0, (y - centre.1) * ahead.1) > 0.0 {
                [(0.0, 0.0); 2]
            } else {
                section.column(off)
            };
            let chord = if px > 0.0 {
                [(0.0, 0.0); 2]
            } else {
                section.column(py)
            };
            span(&arc) + span(&chord) - overlap(&arc, &chord)
        },
    )
}

/// A D-shaped block, a disc up to a chord, its arc's and its chord's top
/// creases rounded where they meet across a convex turn: the arc turns
/// about the drum's axis and the chord is straight, so the bands' sections
/// differ in the plane bisecting the turn, and each runs out past the
/// vertex as it does alone, the cuts taking their union. With the chord off
/// the centre and through it, where the arc's cap lies in the chord's
/// wall, both orders take that union, measured within a box round one
/// corner a column at a time. No junction was taken there: the later band
/// ran on through the earlier, and with the chord through the centre the
/// arc built after the chord left part of the union standing.
#[test]
fn rounds_of_a_d_shaped_block_take_both_bands_as_alone() {
    let r = 0.5;
    let half = 1.2;
    for c in [0.5_f64, 0.0] {
        let mut model = ogeom_topo::Model::new();
        let (solid, at, edges) = d_shaped_block(&mut model, c, 0.0, 0.0);
        let want = d_shaped_corner_removed((at.x, at.y), c, r, half);
        for got in removed_near(&model, &solid, &edges, r, at, half) {
            assert!(
                (got - want).abs() < 5e-6,
                "chord at {c}: {got} against {want}"
            );
        }
    }
}

/// What rounding the arc and the chord of [`d_shaped_block`] takes within
/// the square of half side `half` round its corner `at`, the chord's wall
/// leaning `d` and the drum drafted by `k`, a column at a time: the union
/// of what each band takes as it does alone, in the material. The arc's
/// band is its torus band run on past the vertex to [`arc_stop`]; the
/// chord's is its straight band run out past the vertex
/// across the drum, in the material within the drum's wall at each height.
fn leaning_d_corner_removed(at: (f64, f64), c: f64, (d, k): (f64, f64), r: f64, half: f64) -> f64 {
    let (centre, big) = ((4.0, 4.0), 2.5);
    let arc_band = Section::new(0.0, (1.0, k), r);
    let ahead = arc_stop(c, big, arc_band.centre.0);
    let chord_band = Section::new(0.0, (1.0, d), r);
    integral(
        (at.0 - half, at.0 + half),
        (at.1 - half, at.1 + half),
        2e3,
        |x, y| {
            let off = (x - centre.0).hypot(y - centre.1) - big;
            // The heights the drum holds at this column: those where its
            // radius, `big + k (2 - z)`, reaches `off` past the top's.
            let drum = if k == 0.0 {
                if off > 0.0 { (0.0, 0.0) } else { (0.0, 2.0) }
            } else if k > 0.0 {
                (0.0, (2.0 - off / k).min(2.0))
            } else {
                ((2.0 - off / k).max(0.0), 2.0)
            };
            let py = y - at.1;
            let within = |runs: [(f64, f64); 2]| -> Vec<(f64, f64)> {
                runs.iter()
                    .map(|&(p, q)| (p.max(0.0), q.min(2.0)))
                    .collect()
            };
            let arc = if (x - centre.0).mul_add(ahead.0, (y - centre.1) * ahead.1) > 0.0 {
                Vec::new()
            } else {
                within(arc_band.column(off))
            };
            let chord = meet(&within(chord_band.column(py)), &[drum]);
            let past = [chord_band.past(py, 0.0)];
            span(&arc) + span(&chord) - overlap(&arc, &chord) - overlap(&arc, &past)
        },
    )
}

/// A D-shaped block with its chord through the centre and the chord's wall
/// overhanging, its arc's and its chord's top creases rounded in both
/// orders. The chord's wall cuts the arc's torus band along a section that
/// leaves the wall tangent to its top edge, where the band is tangent to
/// the top, and to its seam with the drum, where the band is tangent to the
/// drum. Named arc first, the arc's cut met each edge anywhere in the
/// valley the section's fit leaves along a tangent, 5.7e-4 and 1.3e-4 off
/// the edges' own touches, the wall's pieces did not close with the band's,
/// and the request refused.
#[test]
fn rounds_of_a_d_shaped_block_with_an_overhanging_chord_close_in_both_orders() {
    let r = 0.5;
    let half = 1.2;
    for d in [-0.1_f64, -0.05] {
        let mut model = ogeom_topo::Model::new();
        let (solid, at, edges) = d_shaped_block(&mut model, 0.0, d, 0.0);
        let want = leaning_d_corner_removed((at.x, at.y), 0.0, (d, 0.0), r, half);
        for got in removed_near(&model, &solid, &edges, r, at, half) {
            assert!(
                (got - want).abs() < 5e-6,
                "chord wall leaning {d}: {got} against {want}"
            );
        }
    }
}

/// A D-shaped block whose drum is drafted (a cone), its foot out or in,
/// its chord's wall square or leaning, its arc's and its chord's top
/// creases rounded in both orders. The arc's round is the torus band of
/// its section, a cone's, and the turn is convex, so each band runs out as
/// it does alone and the cuts take their union, the chord's band run out
/// across the cone. The arc's round was fitted, and the request refused in
/// both orders. With the chord off the centre, its band run out across the
/// cone met the arc's on the same cone spelled on another frame, which the
/// boolean could only image by a fit that missed: named chord first, the
/// junction refused and the bands met as built, 1.1e-4 short with the foot
/// out and the chord leaning, 7e-6 and 2.3e-4 over with the foot in.
#[test]
fn rounds_of_a_d_shaped_block_on_a_drafted_drum_take_both_bands_as_alone() {
    let r = 0.5;
    let half = 1.2;
    for (c, d, k) in [
        (0.5_f64, 0.0_f64, 0.1_f64),
        (0.0, 0.0, 0.1),
        (0.0, 0.0, -0.1),
        (0.5, 0.1, 0.1),
        (0.5, 0.0, -0.1),
        (0.5, 0.1, -0.1),
    ] {
        let mut model = ogeom_topo::Model::new();
        let (solid, at, edges) = d_shaped_block(&mut model, c, d, k);
        let want = leaning_d_corner_removed((at.x, at.y), c, (d, k), r, half);
        for got in removed_near(&model, &solid, &edges, r, at, half) {
            assert!(
                (got - want).abs() < 5e-6,
                "chord at {c} leaning {d}, drum drafted {k}: {got} against {want}"
            );
        }
    }
}

/// Composite 8-point Gauss-Legendre on intervals with no branch changes.
/// Callers split at analytic contact boundaries and transform square-root
/// endpoints before integrating.
fn gauss(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    const X: [f64; 4] = [
        0.183_434_642_495_649_8,
        0.525_532_409_916_329,
        0.796_666_477_413_626_7,
        0.960_289_856_497_536_3,
    ];
    const W: [f64; 4] = [
        0.362_683_783_378_362,
        0.313_706_645_877_887_3,
        0.222_381_034_453_374_5,
        0.101_228_536_290_376_3,
    ];
    let h = (b - a) / 64.0;
    (0..64)
        .map(|k| {
            let mid = a + h * (f64::from(k) + 0.5);
            X.iter()
                .zip(W)
                .map(|(x, w)| w * (f(mid + 0.5 * h * x) + f(mid - 0.5 * h * x)))
                .sum::<f64>()
                * 0.5
                * h
        })
        .sum()
}

fn gauss_split(f: impl Fn(f64) -> f64, mut breaks: Vec<f64>) -> f64 {
    breaks.sort_by(f64::total_cmp);
    breaks.windows(2).map(|ab| gauss(&f, ab[0], ab[1])).sum()
}

/// A disk clipped by the two straight sides of a rounded box's top loop.
/// The drill is clear of the other two sides. At height `h`, both sides
/// have inset `r - sqrt(r²-h²)` from their original positions.
fn drill_area_in_top_round(
    (width, length): (f64, f64),
    top_radius: f64,
    centre: (f64, f64),
    radius: f64,
    h: f64,
) -> f64 {
    let inset = top_radius - (top_radius * top_radius - h * h).sqrt();
    let y0 = inset.max(centre.1 - radius);
    let y1 = (length - inset).min(centre.1 + radius);
    if y1 <= y0 {
        return 0.0;
    }
    let right = width - inset;
    let a = ((y0 - centre.1) / radius).asin();
    let b = ((y1 - centre.1) / radius).asin();
    let mut breaks = vec![a, b];
    let dx = right - centre.0;
    if dx.abs() < radius {
        // A drill chord meets the straight side at these two y values.
        let dy = (radius * radius - dx * dx).sqrt();
        for y in [centre.1 - dy, centre.1 + dy] {
            if (y0..y1).contains(&y) {
                breaks.push(((y - centre.1) / radius).asin());
            }
        }
    }
    gauss_split(
        |theta| {
            let chord = radius * theta.cos();
            ((centre.0 + chord).min(right) - (centre.0 - chord)).max(0.0) * radius * theta.cos()
        },
        breaks,
    )
}

/// The edges of `shape` on the circle of `radius` about `centre`, square to
/// the `z` axis.
fn circle_edges(
    model: &ogeom_topo::Model,
    shape: &ogeom_topo::Shape,
    centre: Point,
    radius: f64,
) -> Vec<ogeom_topo::Shape> {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let data = model.node(e).unwrap().data().as_edge().unwrap();
            let Some(ogeom_topo::EdgeRepr::Curve3d { curve, .. }) = data.curve3d() else {
                return false;
            };
            let Some(ogeom_geom::Curve::Circle(c)) = model.geometry().curve(*curve) else {
                return false;
            };
            let circle = c.circle();
            circle.centre().distance(centre) < 1e-9 && (circle.radius() - radius).abs() < 1e-9
        })
        .collect()
}

/// What a rim's rounding at `r` takes from the solid, the section between
/// the wall, the cap and the tube turned about the axis: its area
/// `(1 - pi/4) r^2`, its centroid `r / (6 (1 - pi/4))` from the tube's
/// centre toward the corner, which stands at `tube`, out from the axis or
/// in toward it.
fn rim_removed(tube: f64, r: f64, outward: bool) -> f64 {
    let area = (1.0 - core::f64::consts::FRAC_PI_4) * r * r;
    let offset = r / (6.0 * (1.0 - core::f64::consts::FRAC_PI_4));
    let at = if outward {
        tube + offset
    } else {
        tube - offset
    };
    core::f64::consts::TAU * at * area
}

/// Rims between coaxial planes and cylinders have the volume obtained by
/// revolving their circular spandrel section, including intersections
/// with cavities and transverse drills.
/// The cases are a cylinder's rims, both on one wall, a bore's through a
/// block, a tube's four (its outer and inner rims each on the one annulus),
/// a wall with a cross hole well below its rim, a cup's rim clear of its
/// cavity, and a drum as tall as the rounding (its wall's ring on the far
/// rim; the closed-form volume holds). Where a hole in the wall comes up
/// into the band the round takes, or a cavity into the box it can take up,
/// the engine hands the request back and the general fuse cuts it. No case
/// is cut by trimming any more.
#[test]
fn rims_are_cut_by_trimming_to_the_general_fuse_s_solid() {
    use core::f64::consts::PI;
    type Made = (ogeom_topo::Shape, Vec<ogeom_topo::Shape>, f64, Option<f64>);
    type Case<'a> = &'a dyn Fn(&mut ogeom_topo::Model) -> Made;
    let cylinder = |model: &mut ogeom_topo::Model, at: Point, axis: Direction, r: f64, h: f64| {
        let frame = Frame::new(
            at,
            axis,
            Direction::new(axis.vector().cross(Vector::new(0.3, 0.5, 0.7)), T).unwrap(),
            T,
        )
        .unwrap();
        ogeom_algo::make_cylinder(model, frame, r, h, T)
            .unwrap()
            .shape
    };
    let drum = |model: &mut ogeom_topo::Model, h: f64, both: bool| -> Made {
        let solid = cylinder(model, Point::ORIGIN, Direction::Z, 5.0, h);
        let mut edges = circle_edges(model, &solid, Point::new(0.0, 0.0, h), 5.0);
        if both {
            edges.extend(circle_edges(model, &solid, Point::ORIGIN, 5.0));
        }
        let rims = if both { 2.0 } else { 1.0 };
        let want = PI * 25.0 * h - rims * rim_removed(4.0, 1.0, true);
        (solid, edges, 1.0, Some(want))
    };
    let bore = |model: &mut ogeom_topo::Model| -> Made {
        let block = ogeom_algo::make_box(model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let drill = cylinder(model, Point::new(10.0, 10.0, -1.0), Direction::Z, 4.0, 12.0);
        let solid = ogeom_bool::cut(model, &block, &drill, T).unwrap().shape;
        let edges = circle_edges(model, &solid, Point::new(10.0, 10.0, 10.0), 4.0);
        let want = 4000.0 - PI * 160.0 - rim_removed(5.0, 1.0, false);
        (solid, edges, 1.0, Some(want))
    };
    let tube = |model: &mut ogeom_topo::Model| -> Made {
        let outer = cylinder(model, Point::ORIGIN, Direction::Z, 5.0, 10.0);
        let inner = cylinder(model, Point::new(0.0, 0.0, -1.0), Direction::Z, 3.0, 12.0);
        let solid = ogeom_bool::cut(model, &outer, &inner, T).unwrap().shape;
        let mut edges = Vec::new();
        for z in [0.0, 10.0] {
            for r in [3.0, 5.0] {
                edges.extend(circle_edges(model, &solid, Point::new(0.0, 0.0, z), r));
            }
        }
        let want = PI * 16.0 * 10.0
            - 2.0 * rim_removed(4.5, 0.5, true)
            - 2.0 * rim_removed(3.5, 0.5, false);
        (solid, edges, 0.5, Some(want))
    };
    let crossed = |model: &mut ogeom_topo::Model, height: f64| -> Made {
        let solid = cylinder(model, Point::ORIGIN, Direction::Z, 5.0, 10.0);
        let drill = cylinder(
            model,
            Point::new(-6.0, 0.0, height),
            Direction::X,
            1.0,
            12.0,
        );
        let solid = ogeom_bool::cut(model, &solid, &drill, T).unwrap().shape;
        let edges = circle_edges(model, &solid, Point::new(0.0, 0.0, 10.0), 5.0);
        // Perpendicular cylinders: integrate the drill's vertical chord and
        // the drum's horizontal chord over their common y coordinates.
        let hole = gauss(
            |theta| {
                let y = theta.sin();
                4.0 * (25.0 - y * y).sqrt() * theta.cos().powi(2)
            },
            -PI / 2.0,
            PI / 2.0,
        );
        // Where the drill reaches the rim, add back the intersection of its
        // void with the spandrel: annuli from the rolling-ball surface to
        // the drum wall, restricted to the drill's circular sections.
        let shared = if height == 8.5 {
            // z = height + cos(theta): the drill chord in y vanishes
            // smoothly at the upper contact, where theta = 0.
            gauss(
                |theta| {
                    let z = height + theta.cos();
                    let y1 = theta.sin();
                    let inner = 4.0 + (1.0 - (z - 9.0).powi(2)).sqrt();
                    gauss(
                        |t| {
                            let y = y1 * t;
                            2.0 * ((25.0 - y * y).sqrt() - (inner * inner - y * y).sqrt()) * y1
                        },
                        -1.0,
                        1.0,
                    ) * theta.sin()
                },
                0.0,
                (9.0 - height).acos(),
            )
        } else {
            0.0
        };
        let want = PI * 250.0 - hole - rim_removed(4.0, 1.0, true) + shared;
        (solid, edges, 1.0, Some(want))
    };
    // A drum hollowed from below, its cavity coming up to `top` at `wall`
    // from the axis.
    let cup = |model: &mut ogeom_topo::Model, wall: f64, top: f64| -> Made {
        let solid = cylinder(model, Point::ORIGIN, Direction::Z, 5.0, 10.0);
        let cavity = cylinder(
            model,
            Point::new(0.0, 0.0, -1.0),
            Direction::Z,
            wall,
            top + 1.0,
        );
        let solid = ogeom_bool::cut(model, &solid, &cavity, T).unwrap().shape;
        let edges = circle_edges(model, &solid, Point::new(0.0, 0.0, 10.0), 5.0);
        // Add back any rim spandrel already removed by the cavity. At z its
        // inner radius is the rolling ball's wall contact curve.
        // The cavity first intersects the spandrel where its radius meets
        // the inner arc, at z = 9 + sqrt(1 - (wall - 4)²).
        let contact = if wall <= 4.0 {
            10.0
        } else if wall >= 5.0 {
            9.0
        } else {
            9.0 + (1.0 - (wall - 4.0).powi(2)).sqrt()
        };
        let shared = if top > contact {
            gauss(
                |theta| {
                    let inner = 4.0 + theta.cos();
                    PI * (wall * wall - inner * inner) * theta.cos()
                },
                (contact - 9.0).asin(),
                (top.min(10.0) - 9.0).asin(),
            )
        } else {
            0.0
        };
        let want = PI * (250.0 - wall * wall * top) - rim_removed(4.0, 1.0, true) + shared;
        (solid, edges, 1.0, Some(want))
    };
    let cases: [(&str, Case, bool); 9] = [
        ("a drum's rim", &|m| drum(m, 10.0, false), true),
        ("both rims of a short drum", &|m| drum(m, 2.5, true), true),
        ("a bore's rim", &bore, true),
        ("a tube's rims", &tube, true),
        ("a wall with a cross hole below", &|m| crossed(m, 3.0), true),
        (
            "a drum as tall as the rounding",
            &|m| drum(m, 1.0, false),
            true,
        ),
        ("a cross hole into the band", &|m| crossed(m, 8.5), false),
        (
            "a cup's rim clear of its cavity",
            &|m| cup(m, 3.5, 8.0),
            true,
        ),
        (
            "a cup's rim breaking into its cavity",
            &|m| cup(m, 4.7, 9.7),
            false,
        ),
    ];
    for (name, case, by_engine) in cases {
        let mut model = ogeom_topo::Model::new();
        let (solid, edges, radius, closed_form) = case(&mut model);
        let (seated, work) = ogeom_bool::work::counted(|| {
            ogeom_fillet::fillet_edges(&mut model, &solid, &edges, radius, T).unwrap()
        });
        assert_eq!(
            work.face_pairs == 0,
            by_engine,
            "{name}: built with no boolean"
        );
        assert_eq!(work.seated_cuts, 0, "{name}: no cut by trimming");
        assert!(
            ogeom_algo::check(&model, &seated.shape, T)
                .unwrap()
                .is_valid(),
            "{name}: the rounding is valid"
        );
        let got = volume(&model, &seated.shape, 1e-3);
        let exact = closed_form.expect("every rim has an independent volume");
        // The 1e-3 chord's curved-face meshing error is bounded relatively;
        // the 64-panel quadratures are much more accurate than this tolerance.
        assert!(
            (got - exact).abs() <= exact * 1e-7,
            "{name}: volume {got} against the independent value {exact}"
        );
    }
}

/// Concave rims where a plane meets a cylinder add the volume obtained by
/// revolving the spandrel section, accounting for other joined features.
/// The cases are a boss's base on a plate, a blind hole's floor, a stepped shaft's shoulder, two bosses'
/// bases in one union, a bored boss's base, a shelf over the ring's top, a
/// boss whose ring just meets the plate's side and a boss as tall as the
/// rounding (their closed-form volumes hold). Where a shelf hangs into the
/// box the round can take up, the engine hands the request back and the
/// general fuse fills it. No case is filled by trimming any more.
#[test]
fn concave_rims_are_filled_by_trimming_to_the_general_fuse_s_solid() {
    use core::f64::consts::PI;
    type Made = (ogeom_topo::Shape, Vec<ogeom_topo::Shape>, Option<f64>);
    type Case<'a> = &'a dyn Fn(&mut ogeom_topo::Model) -> Made;
    let cylinder = |model: &mut ogeom_topo::Model, at: Point, r: f64, h: f64| {
        let frame = Frame::new(at, Direction::Z, Direction::X, T).unwrap();
        ogeom_algo::make_cylinder(model, frame, r, h, T)
            .unwrap()
            .shape
    };
    let plate = |model: &mut ogeom_topo::Model, x: f64| {
        ogeom_algo::make_box(model, Frame::WORLD, (x, 20.0, 5.0), T)
            .unwrap()
            .shape
    };
    // A boss of radius 3 and height `h` at `(x, 10)` on a 20 x 20 plate 5
    // thick, its base rounded at 1.
    let boss = |model: &mut ogeom_topo::Model, x: f64, h: f64| -> Made {
        let base = plate(model, 20.0);
        let post = cylinder(model, Point::new(x, 10.0, 5.0), 3.0, h);
        let solid = ogeom_bool::fuse(model, &base, &post, T).unwrap().shape;
        let edges = circle_edges(model, &solid, Point::new(x, 10.0, 5.0), 3.0);
        let want = 2000.0 + PI * 9.0 * h + rim_removed(4.0, 1.0, false);
        (solid, edges, Some(want))
    };
    let blind = |model: &mut ogeom_topo::Model| -> Made {
        let block = ogeom_algo::make_box(model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let drill = cylinder(model, Point::new(10.0, 10.0, 5.0), 3.0, 6.0);
        let solid = ogeom_bool::cut(model, &block, &drill, T).unwrap().shape;
        let edges = circle_edges(model, &solid, Point::new(10.0, 10.0, 5.0), 3.0);
        let want = 4000.0 - PI * 45.0 + rim_removed(2.0, 1.0, true);
        (solid, edges, Some(want))
    };
    let shaft = |model: &mut ogeom_topo::Model| -> Made {
        let wide = cylinder(model, Point::ORIGIN, 5.0, 5.0);
        let narrow = cylinder(model, Point::new(0.0, 0.0, 5.0), 3.0, 5.0);
        let solid = ogeom_bool::fuse(model, &wide, &narrow, T).unwrap().shape;
        let edges = circle_edges(model, &solid, Point::new(0.0, 0.0, 5.0), 3.0);
        let want = PI * 125.0 + PI * 45.0 + rim_removed(4.0, 1.0, false);
        (solid, edges, Some(want))
    };
    let bosses = |model: &mut ogeom_topo::Model| -> Made {
        let mut solid = plate(model, 30.0);
        let mut edges = Vec::new();
        for x in [8.0, 22.0] {
            let post = cylinder(model, Point::new(x, 10.0, 5.0), 3.0, 5.0);
            solid = ogeom_bool::fuse(model, &solid, &post, T).unwrap().shape;
        }
        for x in [8.0, 22.0] {
            edges.extend(circle_edges(model, &solid, Point::new(x, 10.0, 5.0), 3.0));
        }
        let want = 3000.0 + 2.0 * (PI * 45.0 + rim_removed(4.0, 1.0, false));
        (solid, edges, Some(want))
    };
    let bored = |model: &mut ogeom_topo::Model| -> Made {
        let (solid, _, _) = boss(model, 10.0, 5.0);
        let drill = cylinder(model, Point::new(10.0, 10.0, -1.0), 1.5, 12.0);
        let solid = ogeom_bool::cut(model, &solid, &drill, T).unwrap().shape;
        let edges = circle_edges(model, &solid, Point::new(10.0, 10.0, 5.0), 3.0);
        let want = 2000.0 + PI * 45.0 - PI * 2.25 * 10.0 + rim_removed(4.0, 1.0, false);
        (solid, edges, Some(want))
    };
    // A boss with a shelf on a post beside it, the shelf's underside `at`
    // over the plate.
    let shelved = |model: &mut ogeom_topo::Model, at: f64| -> Made {
        let (mut solid, _, _) = boss(model, 10.0, 5.0);
        for (corner, size) in [
            (Point::new(16.0, 8.0, 5.0), (2.0, 4.0, 3.0)),
            (Point::new(13.5, 8.0, 5.0 + at), (4.5, 4.0, 1.0)),
        ] {
            let frame = Frame::new(corner, Direction::Z, Direction::X, T).unwrap();
            let part = ogeom_algo::make_box(model, frame, size, T).unwrap().shape;
            solid = ogeom_bool::fuse(model, &solid, &part, T).unwrap().shape;
        }
        let edges = circle_edges(model, &solid, Point::new(10.0, 10.0, 5.0), 3.0);
        // The shelf overlaps the supporting post by 2 x 4 x (8 - (5 + at))
        // when it sits below the post's top. It is clear of the boss.
        let overlap = 8.0 * (3.0 - at).clamp(0.0, 1.0);
        let want = 2000.0 + PI * 45.0 + 42.0 - overlap + rim_removed(4.0, 1.0, false);
        (solid, edges, Some(want))
    };
    let cases: [(&str, Case, bool); 9] = [
        ("a boss's base", &|m| boss(m, 10.0, 5.0), true),
        ("a blind hole's floor", &blind, true),
        ("a shaft's shoulder", &shaft, true),
        ("two bosses' bases", &bosses, true),
        ("a bored boss's base", &bored, true),
        ("a shelf over the ring's top", &|m| shelved(m, 2.0), true),
        (
            "a boss whose ring meets the plate's side",
            &|m| boss(m, 4.0, 5.0),
            true,
        ),
        ("a shelf into the ring", &|m| shelved(m, 0.5), false),
        (
            "a boss as tall as the rounding",
            &|m| boss(m, 10.0, 1.0),
            true,
        ),
    ];
    for (name, case, by_engine) in cases {
        let mut model = ogeom_topo::Model::new();
        let (solid, edges, closed_form) = case(&mut model);
        let (seated, work) = ogeom_bool::work::counted(|| {
            ogeom_fillet::fillet_edges(&mut model, &solid, &edges, 1.0, T).unwrap()
        });
        assert_eq!(
            work.face_pairs == 0,
            by_engine,
            "{name}: built with no boolean"
        );
        assert_eq!(work.seated_fills, 0, "{name}: no fill by trimming");
        assert!(
            ogeom_algo::check(&model, &seated.shape, T)
                .unwrap()
                .is_valid(),
            "{name}: the rounding is valid"
        );
        let got = volume(&model, &seated.shape, 1e-3);
        let exact = closed_form.expect("every rim has an independent volume");
        // The 1e-3 chord's curved-face meshing error is bounded relatively;
        // the closed forms and Pappus sections are evaluated in f64.
        assert!(
            (got - exact).abs() <= exact * 1e-7,
            "{name}: volume {got} against the independent value {exact}"
        );
    }
}

/// Non-overlapping rims and straight edges take their independently
/// calculated spandrel volumes in one request. A boss's top rim with the
/// plate's top loop, with its bottom loop or with one bottom edge, and a
/// tube boss's two top rims with the plate's bottom loop, go through the
/// engine's trimmed assembly whole, with no boolean. For a union, a boss's
/// base fused on with a rib's concave base edges, the rim's wedge is
/// trimmed out first and the rest cut from what is left, the trim keeping
/// only the faces that cut sets aside.
#[test]
fn rims_grouped_with_other_wedges_are_trimmed_before_the_rest_is_cut() {
    use core::f64::consts::{FRAC_PI_4, PI};
    type Made = (ogeom_topo::Shape, Vec<ogeom_topo::Shape>, f64, f64);
    let section = 1.0 - FRAC_PI_4;
    let cell = 2.0 * (5.0 / 6.0 - FRAC_PI_4);
    type Case<'a> = &'a dyn Fn(&mut ogeom_topo::Model) -> Made;
    // A 20 x 20 plate 5 thick with a boss of radius `r` and height 5 on it,
    // bored through at radius `bore` where one is given, the seams off the
    // plate's axes.
    let seam = Direction::new(Vector::new(-0.5, 0.3, 0.0), T).unwrap();
    let boss = |model: &mut ogeom_topo::Model, r: f64, bore: Option<f64>| {
        let plate = ogeom_algo::make_box(model, Frame::WORLD, (20.0, 20.0, 5.0), T)
            .unwrap()
            .shape;
        let frame = Frame::new(Point::new(10.0, 10.0, 5.0), Direction::Z, seam, T).unwrap();
        let post = ogeom_algo::make_cylinder(model, frame, r, 5.0, T)
            .unwrap()
            .shape;
        let mut solid = ogeom_bool::fuse(model, &plate, &post, T).unwrap().shape;
        if let Some(bore) = bore {
            let frame = Frame::new(Point::new(10.0, 10.0, -1.0), Direction::Z, seam, T).unwrap();
            let drill = ogeom_algo::make_cylinder(model, frame, bore, 12.0, T)
                .unwrap()
                .shape;
            solid = ogeom_bool::cut(model, &solid, &drill, T).unwrap().shape;
        }
        solid
    };
    // The plate's sides at height `z`: its edges between two of its corners.
    let loop_at = |model: &ogeom_topo::Model, solid: &ogeom_topo::Shape, z: f64| {
        let corner = |v: &ogeom_topo::Shape| {
            let p = model.node(v).unwrap().data().as_vertex().unwrap().point;
            [p.x, p.y]
                .iter()
                .all(|c| c.abs() < 1e-9 || (c - 20.0).abs() < 1e-9)
                && (p.z - z).abs() < 1e-9
        };
        explore_unique(model, solid, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                ogeom_algo::edge_vertices(model, e)
                    .unwrap()
                    .is_some_and(|(a, b)| corner(&a) && corner(&b))
            })
            .collect::<Vec<_>>()
    };
    let top = |model: &mut ogeom_topo::Model, z: f64, all: bool| -> Made {
        let solid = boss(model, 3.0, None);
        let mut edges = circle_edges(model, &solid, Point::new(10.0, 10.0, 10.0), 3.0);
        let straight = loop_at(model, &solid, z);
        let take = if all { straight.len() } else { 1 };
        edges.extend(straight.into_iter().take(take));
        let plate = if all {
            80.0 * section - 4.0 * cell
        } else {
            20.0 * section
        };
        let want = 2000.0 + PI * 45.0 - rim_removed(2.0, 1.0, true) - plate;
        (solid, edges, 1.0, want)
    };
    let tube = |model: &mut ogeom_topo::Model| -> Made {
        let solid = boss(model, 4.0, Some(2.0));
        let mut edges = Vec::new();
        for r in [4.0, 2.0] {
            edges.extend(circle_edges(model, &solid, Point::new(10.0, 10.0, 10.0), r));
        }
        edges.extend(loop_at(model, &solid, 0.0));
        let s = section * 0.25;
        let c = cell * 0.125;
        let want = 2000.0 + PI * 40.0
            - rim_removed(3.5, 0.5, true)
            - rim_removed(2.5, 0.5, false)
            - (80.0 * s - 4.0 * c);
        (solid, edges, 0.5, want)
    };
    // A boss of radius 3 at (8, 10) on a 30 x 20 plate and a rib 2 wide
    // standing on it from x = 20, the boss's base and the rib's long base
    // edges rounded: concave, fused on in one union.
    let ribbed = |model: &mut ogeom_topo::Model| -> Made {
        let plate = ogeom_algo::make_box(model, Frame::WORLD, (30.0, 20.0, 5.0), T)
            .unwrap()
            .shape;
        let frame = Frame::new(Point::new(8.0, 10.0, 5.0), Direction::Z, seam, T).unwrap();
        let post = ogeom_algo::make_cylinder(model, frame, 3.0, 5.0, T)
            .unwrap()
            .shape;
        let frame = Frame::new(Point::new(20.0, 2.0, 5.0), Direction::Z, Direction::X, T).unwrap();
        let rib = ogeom_algo::make_box(model, frame, (2.0, 16.0, 5.0), T)
            .unwrap()
            .shape;
        let solid = ogeom_bool::fuse(model, &plate, &post, T).unwrap().shape;
        let solid = ogeom_bool::fuse(model, &solid, &rib, T).unwrap().shape;
        let mut edges = circle_edges(model, &solid, Point::new(8.0, 10.0, 5.0), 3.0);
        for x in [20.0, 22.0] {
            edges.extend(
                explore_unique(model, &solid, ShapeType::Edge)
                    .unwrap()
                    .into_iter()
                    .filter(|e| {
                        ogeom_algo::edge_vertices(model, e)
                            .unwrap()
                            .is_some_and(|(a, b)| {
                                [a, b].iter().all(|v| {
                                    let p =
                                        model.node(v).unwrap().data().as_vertex().unwrap().point;
                                    (p.x - x).abs() < 1e-9 && (p.z - 5.0).abs() < 1e-9
                                })
                            })
                    }),
            );
        }
        let want = 3000.0 + PI * 45.0 + 160.0 + rim_removed(4.0, 1.0, false) + 32.0 * section;
        (solid, edges, 1.0, want)
    };
    let cases: [(&str, Case, bool); 5] = [
        (
            "a boss's top rim and the plate's top loop",
            &|m| top(m, 5.0, true),
            false,
        ),
        (
            "a boss's top rim and the plate's bottom loop",
            &|m| top(m, 0.0, true),
            false,
        ),
        (
            "a boss's top rim and one bottom edge",
            &|m| top(m, 0.0, false),
            false,
        ),
        (
            "a tube boss's top rims and the plate's bottom loop",
            &tube,
            false,
        ),
        ("a boss's base and a rib's base edges", &ribbed, true),
    ];
    for (name, case, fill) in cases {
        let mut model = ogeom_topo::Model::new();
        let (solid, edges, radius, want) = case(&mut model);
        let (seated, work) = ogeom_bool::work::counted(|| {
            ogeom_fillet::fillet_edges(&mut model, &solid, &edges, radius, T).unwrap()
        });
        let trimmed = if fill {
            work.seated_fills
        } else {
            work.seated_cuts
        };
        // One bottom edge beside the boss's top rim: both rounds go through
        // the engine's trimmed assembly, with no boolean.
        if !fill {
            assert_eq!(work.face_pairs, 0, "{name}: built with no boolean");
            assert_eq!(trimmed, 0, "{name}: nothing trimmed");
        } else {
            assert!(trimmed > 0, "{name}: the rims are trimmed");
        }
        assert!(
            ogeom_algo::check(&model, &seated.shape, T)
                .unwrap()
                .is_valid(),
            "{name}: the rounding is valid"
        );
        let got = volume(&model, &seated.shape, 1e-3);
        // Pappus and the prism corner cell are exact; 1e-7 relative covers
        // the 1e-3 chord used to mesh the cylindrical and toroidal faces.
        assert!(
            (got - want).abs() <= want * 1e-7,
            "{name}: volume {got} against the independent value {want}"
        );
    }
}

/// Planar bands, including ones ending on a corner or a hole, have the
/// volume of their section swept along each edge with corner overlaps
/// accounted for. Where a band's rail leaves a vertex along another edge closely enough
/// that the general fuse splits the two where they part (a corner rounded
/// at a hundredth, a band ending on an edge stating a millimetre's
/// thousandth), the cut is left to the general fuse. A solid the blend
/// engine takes (the cube, the corners, the bands across, a slab rounded at
/// half its thickness, its sides consumed, all of planes, and the holed
/// block, whose bands run along planes beside a bore and go through its
/// trimmed assembly) is cut by neither: it is built in one pass, the same
/// solid either way.
#[test]
fn seated_bands_are_cut_by_trimming_to_the_general_fuse_s_solid() {
    use core::f64::consts::{FRAC_PI_4, PI};
    type Made = (ogeom_topo::Shape, Vec<ogeom_topo::Shape>, f64, f64);
    let section = 1.0 - FRAC_PI_4;
    // Steiner's formula for an inner box swept by a ball: inner volume,
    // planar slabs, quarter-cylinders along its edges, and eight ball octants.
    let rounded_box = |size: (f64, f64, f64), r: f64| {
        let (x, y, z) = (size.0 - 2.0 * r, size.1 - 2.0 * r, size.2 - 2.0 * r);
        x * y * z
            + 2.0 * (x * y + x * z + y * z) * r
            + PI * r * r * (x + y + z)
            + 4.0 * PI * r * r * r / 3.0
    };
    type Case<'a> = &'a dyn Fn(&mut ogeom_topo::Model) -> Made;
    let cube = |model: &mut ogeom_topo::Model| -> Made {
        let solid = ogeom_algo::make_box(model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let edges = explore_unique(model, &solid, ShapeType::Edge).unwrap();
        let want = rounded_box((10.0, 10.0, 10.0), 1.0);
        (solid, edges, 1.0, want)
    };
    let corner = |model: &mut ogeom_topo::Model, radius: f64| -> Made {
        let solid = ogeom_algo::make_box(model, Frame::WORLD, (2.0, 2.0, 1.0), T)
            .unwrap()
            .shape;
        let edges = [
            Point::new(2.0, 1.0, 1.0),
            Point::new(2.0, 2.0, 0.5),
            Point::new(1.0, 2.0, 1.0),
        ]
        .map(|at| edge_through(model, &solid, at))
        .to_vec();
        // The three mutually perpendicular rounds meet in a spherical
        // octant. The corner correction is the cubic term of Steiner's
        // rounded-box volume divided among its eight identical corners.
        let cell = radius.powi(3) * (2.0 - 7.0 * PI / 12.0);
        let want = 4.0 - 5.0 * section * radius * radius + cell;
        (solid, edges, radius, want)
    };
    let across = |model: &mut ogeom_topo::Model, loose: bool| -> Made {
        let solid = ogeom_algo::make_box(model, Frame::WORLD, (10.0, 6.0, 4.0), T)
            .unwrap()
            .shape;
        if loose {
            // The upright edge one band ends on states a thousandth.
            let upright = edge_through(model, &solid, Point::new(10.0, 0.0, 2.0));
            model
                .widen(&upright, ogeom_core::Tolerance::new(1e-3).unwrap())
                .unwrap();
        }
        let edges = [Point::new(5.0, 0.0, 4.0), Point::new(5.0, 6.0, 0.0)]
            .map(|at| edge_through(model, &solid, at))
            .to_vec();
        (solid, edges, 1.0, 240.0 - 20.0 * section)
    };
    // A 10 x 6 plate 4 thick through which a hole of radius 1 runs, its
    // rims each one closed edge: the top the bands sit on walks one.
    let holed = |model: &mut ogeom_topo::Model| -> Made {
        let outline = ogeom_algo::make_polygon(
            model,
            &[
                Point::new(0.0, 0.0, 0.0),
                Point::new(10.0, 0.0, 0.0),
                Point::new(10.0, 6.0, 0.0),
                Point::new(0.0, 6.0, 0.0),
            ],
            true,
            T,
        )
        .unwrap()
        .shape;
        let rim = Circle::new(
            Frame::new(Point::new(5.0, 3.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
            1.0,
            T,
        )
        .unwrap();
        let rim = ogeom_algo::make_edge(
            model,
            ogeom_geom::Curve::Circle(CircleCurve::new(rim)),
            (0.0, core::f64::consts::TAU),
            T,
        )
        .unwrap()
        .shape;
        let hole = ogeom_algo::make_wire(model, &[rim], T).unwrap().shape;
        let plane = ogeom_geom::PlaneSurface::new(Plane::new(Frame::WORLD));
        let face = ogeom_algo::make_face(model, plane.into(), &[outline, hole], T)
            .unwrap()
            .shape;
        let prism = ogeom_algo::make_prism(model, &face, Vector::new(0.0, 0.0, 4.0), T)
            .unwrap()
            .shape;
        // Stated in place, its top no copy of its bottom moved up.
        let solid = ogeom_algo::baked_shape(model, &prism, T).unwrap().shape;
        let edges = [Point::new(5.0, 0.0, 4.0), Point::new(5.0, 6.0, 4.0)]
            .map(|at| edge_through(model, &solid, at))
            .to_vec();
        (solid, edges, 1.0, 240.0 - 4.0 * PI - 20.0 * section)
    };
    let slab = |model: &mut ogeom_topo::Model| -> Made {
        let solid = ogeom_algo::make_box(model, Frame::WORLD, (3.0, 2.0, 1.0), T)
            .unwrap()
            .shape;
        let edges = explore_unique(model, &solid, ShapeType::Edge).unwrap();
        let want = rounded_box((3.0, 2.0, 1.0), 0.5);
        (solid, edges, 0.5, want)
    };
    let cases: [(&str, Case, bool); 7] = [
        ("cube", &cube, false),
        ("corner at a tenth", &|m| corner(m, 0.1), false),
        ("corner at a hundredth", &|m| corner(m, 0.01), false),
        ("across", &|m| across(m, false), false),
        ("across a loose edge", &|m| across(m, true), false),
        ("slab", &slab, false),
        ("holed", &holed, false),
    ];
    for (name, case, trimmed) in cases {
        let mut model = ogeom_topo::Model::new();
        let (solid, edges, radius, want) = case(&mut model);
        let (seated, work) = ogeom_bool::work::counted(|| {
            ogeom_fillet::fillet_edges(&mut model, &solid, &edges, radius, T).unwrap()
        });
        assert_eq!(work.seated_cuts > 0, trimmed, "{name}: cut by trimming");
        assert!(
            ogeom_algo::check(&model, &seated.shape, T)
                .unwrap()
                .is_valid(),
            "{name}: the rounding is valid"
        );
        let got = volume(&model, &seated.shape, 1e-3);
        // Steiner's formula and the straight spandrel are exact; 1e-7
        // relative covers tessellation of the rounded faces at chord 1e-3.
        assert!(
            (got - want).abs() <= want * 1e-7,
            "{name}: volume {got} against the independent value {want}"
        );
    }
}

/// A block whose top `z <= 2` has an outside corner at `(0, 0, 2)` of
/// `plan` radians, between a crease along `x` and one along `(cos plan,
/// sin plan)`, the material between them. Each wall leans by its `leans`
/// entry, its outward normal `(n, k)` normalised for the outward unit `n`
/// square to its crease in plan, so a positive lean puts its foot out.
/// Square walls end each crease six along. Returns the block and the two creases,
/// the one along `x` first.
fn outside_corner(
    model: &mut ogeom_topo::Model,
    plan: f64,
    leans: [f64; 2],
) -> (ogeom_topo::Shape, [ogeom_topo::Shape; 2]) {
    let (s, c) = plan.sin_cos();
    // The plan a unit out past both creases, cut back to the walls.
    let outline = [
        Point::new(-(1.0 + c) / s, -1.0, 0.0),
        Point::new(6.0, -1.0, 0.0),
        Point::new(6.0, 6.0f64.mul_add(-c, 6.0) / s, 0.0),
        Point::new(6.0f64.mul_add(c, -s), 6.0f64.mul_add(s, c), 0.0),
    ];
    let mut block = polygon_prism(model, &outline, Direction::Z, Vector::new(0.0, 0.0, 2.0));
    let top = Point::new(0.0, 0.0, 2.0);
    let along = Vector::new(c, s, 0.0);
    for normal in [
        Vector::new(0.0, -1.0, leans[0]),
        Vector::new(-s, c, leans[1]),
    ] {
        let normal = normal / normal.magnitude();
        let past = past_plane(model, top, normal);
        let cut = ogeom_bool::cut(model, &block, &past, T).unwrap().shape;
        block = explore_unique(model, &cut, ShapeType::Solid).unwrap()[0].clone();
    }
    let edges = [
        placed_edge_through(model, &block, Point::new(3.0, 0.0, 2.0)),
        placed_edge_through(model, &block, top + along * 3.0),
    ];
    (block, edges)
}

/// The runs of `u` the section of `band` holds at height `z`: within its
/// corner's half-planes, less the ball's chord there.
fn slice(band: &Section, z: f64) -> [(f64, f64); 2] {
    let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
    for &((au, az), b) in &band.sides {
        let room = az.mul_add(-z, b);
        if au.abs() < 1e-15 {
            if room < 0.0 {
                return [(0.0, 0.0); 2];
            }
        } else if au > 0.0 {
            hi = hi.min(room / au);
        } else {
            lo = lo.max(room / au);
        }
    }
    if hi <= lo {
        return [(0.0, 0.0); 2];
    }
    let reach = band
        .radius
        .mul_add(band.radius, -(z - band.centre.1).powi(2));
    if reach <= 0.0 {
        return [(lo, hi), (0.0, 0.0)];
    }
    let h = reach.sqrt();
    [
        (lo, hi.min(band.centre.0 - h)),
        (lo.max(band.centre.0 + h), hi),
    ]
}

/// What rounding both creases of [`outside_corner`] at `r` takes
/// under the union rule: each band alone runs out past the vertex across
/// the other's wall, and the two cuts are combined. Taken a height at a
/// time. There each band alone is its section's runs swept along its
/// crease, each line from where it enters the solid through the other's
/// wall to the square wall six along, a length linear in the line's
/// offset, so each run integrates in closed form. What both take is the
/// product of the two runs' widths, `1 / sin plan` the area of a cell in
/// the two offsets.
fn outside_corner_removed(plan: f64, leans: [f64; 2], r: f64) -> f64 {
    let (s, c) = plan.sin_cos();
    let sections = leans.map(|k| Section::new(0.0, (1.0, k), r));
    let steps = 1_000_000;
    let dz = r * 1.5 / f64::from(steps);
    (0..steps)
        .map(|i| {
            let z = 2.0 - (f64::from(i) + 0.5) * dz;
            let runs = [slice(&sections[0], z), slice(&sections[1], z)];
            let width = |runs: &[(f64, f64); 2]| span(runs);
            let alone = |own: &[(f64, f64); 2], other: f64| -> f64 {
                own.iter()
                    .filter(|(u0, u1)| u1 > u0)
                    .map(|&(u0, u1)| {
                        // The line at offset u enters the solid where
                        // `s t >= -u c - other (2 - z)`.
                        let length = |u: f64| 6.0 - other.mul_add(-(2.0 - z), -u * c) / s;
                        (u1 - u0) * f64::midpoint(length(u0), length(u1))
                    })
                    .sum()
            };
            alone(&runs[0], leans[1]) + alone(&runs[1], leans[0])
                - width(&runs[0]) * width(&runs[1]) / s
        })
        .sum::<f64>()
        * dz
}

/// The two top creases at an outside corner, acute, square or obtuse,
/// rounded in either order, the walls square or leaning alike or apart.
/// Each band runs out past the vertex as it does alone and the cuts are
/// combined, as across any convex turn of two straight creases. With walls
/// leaning alike the bands' end sections reflect into each other across
/// the bisecting plane, and with their feet out each section reaches past
/// that plane under the other's: meeting there left a ridge either band
/// alone takes. Measured against the union a height at a time.
#[test]
fn rounds_at_an_outside_corner_take_both_bands_run_out() {
    let r = 0.5;
    for (degrees, leans) in [
        (60.0_f64, [0.0, 0.0]),
        (60.0, [-0.1, -0.1]),
        (60.0, [0.1, 0.1]),
        (45.0, [0.1, 0.1]),
        (90.0, [0.1, 0.1]),
        (120.0, [0.1, 0.1]),
        (60.0, [0.1, -0.05]),
    ] {
        let plan = degrees.to_radians();
        let want = outside_corner_removed(plan, leans, r);
        let mut model = ogeom_topo::Model::new();
        let (block, edges) = outside_corner(&mut model, plan, leans);
        for got in removed_both_ways(&model, &block, &edges, r) {
            assert!(
                (got - want).abs() < 3e-6,
                "{degrees} degrees, leans {leans:?}: {got} against {want}"
            );
        }
    }
}
/// Bands that run, at a block's corner, into a band an earlier group of the
/// same rounding has cut are trimmed, each cut back at the mitre where the
/// two blends cross, to the solid and history the general fuse gives: a
/// block's top loop, two of its top edges, and the top loop of a block
/// standing in a left-handed frame. Each band takes `(1 - pi/4) r^2` along
/// its edge and gives back, at each corner it shares, the cell past the
/// mitre, `2 r^3 (5/6 - pi/4)`. Where a bore comes up beside a corner into
/// the region a band's wedge lies in, or where a block is too narrow for
/// the mitres at its two ends to stand apart, the cut is left to the
/// general fuse.
#[test]
fn bands_running_into_a_rounded_corner_are_trimmed_to_the_general_fuse_s_solid() {
    use core::f64::consts::{FRAC_PI_4, PI};
    let section = 1.0 - FRAC_PI_4;
    let cell = 2.0 * (5.0 / 6.0 - FRAC_PI_4);
    // The edges of a solid's top, `height` up in `frame`, whose middles
    // `pick` keeps, read in the frame.
    let top = |model: &ogeom_topo::Model,
               solid: &ogeom_topo::Shape,
               frame: Frame,
               height: f64,
               pick: &dyn Fn(Point) -> bool| {
        explore_unique(model, solid, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                let data = model.node(e).unwrap().data().as_edge().unwrap();
                let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d()
                else {
                    unreachable!("an edge of a solid made here has a curve");
                };
                let curve = model.geometry().curve(*curve).unwrap();
                let [p, m, q] = [range.0, f64::midpoint(range.0, range.1), range.1]
                    .map(|t| frame.to_local(ogeom_geom::Curve3d::point_at(curve, t, T).unwrap()));
                [p, m, q].iter().all(|a| (a.z - height).abs() < 1e-9) && pick(m)
            })
            .collect::<Vec<_>>()
    };
    let block = |model: &mut ogeom_topo::Model, frame: Frame, size: (f64, f64, f64)| {
        ogeom_algo::make_box(model, frame, size, T).unwrap().shape
    };
    let all = |_: Point| true;
    type Made = (ogeom_topo::Shape, Vec<ogeom_topo::Shape>, f64, Option<f64>);
    type Case<'a> = &'a dyn Fn(&mut ogeom_topo::Model) -> Made;
    let top_loop = |m: &mut ogeom_topo::Model| -> Made {
        let solid = block(m, Frame::WORLD, (20.0, 20.0, 10.0));
        let edges = top(m, &solid, Frame::WORLD, 10.0, &all);
        (
            solid,
            edges,
            1.0,
            Some(4000.0 - 80.0 * section + 4.0 * cell),
        )
    };
    let two_edges = |m: &mut ogeom_topo::Model| -> Made {
        let pick = |p: Point| p.x.abs() < 1e-9 || p.y.abs() < 1e-9;
        let solid = block(m, Frame::WORLD, (20.0, 20.0, 10.0));
        let edges = top(m, &solid, Frame::WORLD, 10.0, &pick);
        (solid, edges, 1.0, Some(4000.0 - 40.0 * section + cell))
    };
    let left_handed = |m: &mut ogeom_topo::Model| -> Made {
        let frame =
            Frame::from_axes(Point::ORIGIN, Direction::Y, Direction::X, Direction::Z, T).unwrap();
        let solid = block(m, frame, (20.0, 14.0, 10.0));
        let edges = top(m, &solid, frame, 10.0, &all);
        (
            solid,
            edges,
            1.0,
            Some(2800.0 - 68.0 * section + 4.0 * cell),
        )
    };
    let bored = |m: &mut ogeom_topo::Model| -> Made {
        let pick = |p: Point| p.x < 1e-9 || p.y < 1e-9 || p.x > 19.9 || p.y > 19.9;
        let solid = block(m, Frame::WORLD, (20.0, 20.0, 10.0));
        let at = Frame::new(Point::new(19.2, 0.6, -1.0), Direction::Z, Direction::X, T).unwrap();
        let bore = ogeom_algo::make_cylinder(m, at, 0.3, 12.0, T)
            .unwrap()
            .shape;
        let solid = ogeom_bool::cut(m, &solid, &bore, T).unwrap().shape;
        let edges = top(m, &solid, Frame::WORLD, 10.0, &pick);
        let rounded = 4000.0 - 80.0 * section + 4.0 * cell;
        // Below the top round the drill removes its full circular section.
        // Above z = 9, integrate only the part inside the rounded profile.
        let (cx, cy, rd): (f64, f64, f64) = (19.2, 0.6, 0.3);
        let dx = 20.0 - cx;
        // The inset d = 1 - cos(theta) reaches the drill's near y edge,
        // the right side's first/last chord contacts and its far y edge.
        // At d = (dx + cy ± sqrt(2rd² - (dx - cy)²))/2 the right-side
        // contact meets the moving lower y boundary.
        let crossing = (2.0 * rd * rd - (dx - cy).powi(2)).sqrt();
        let mut breaks = vec![0.0, core::f64::consts::FRAC_PI_2];
        for d in [
            cy - rd,
            dx - rd,
            (dx + cy - crossing) / 2.0,
            dx,
            (dx + cy + crossing) / 2.0,
            cy + rd,
        ] {
            breaks.push((1.0 - d).acos());
        }
        let drill = 9.0 * PI * rd * rd
            + gauss_split(
                |theta| {
                    drill_area_in_top_round((20.0, 20.0), 1.0, (cx, cy), rd, theta.sin())
                        * theta.cos()
                },
                breaks,
            );
        (solid, edges, 1.0, Some(rounded - drill))
    };
    let narrow = |m: &mut ogeom_topo::Model| -> Made {
        let solid = block(m, Frame::WORLD, (3.0, 20.0, 10.0));
        let edges = top(m, &solid, Frame::WORLD, 10.0, &all);
        // At h above the lower tangent plane, opposite bands meet at a
        // mitre. Their section is 2q(17 + 2q), q = sqrt(r²-h²).
        // The integrals of q and q² from 0 to r are πr²/4 and 2r³/3.
        let r = 1.5;
        let cap = 8.5 * PI * r * r + 8.0 * r.powi(3) / 3.0;
        (solid, edges, r, Some(3.0 * 20.0 * (10.0 - r) + cap))
    };
    let cases: [(&str, Case, bool); 5] = [
        ("a block's top loop", &top_loop, true),
        ("two of a block's top edges", &two_edges, true),
        ("a left-handed block's top loop", &left_handed, true),
        ("a bore beside a corner", &bored, false),
        ("a block too narrow for its mitres", &narrow, false),
    ];
    for (name, case, trimmed) in cases {
        let mut model = ogeom_topo::Model::new();
        let (solid, edges, radius, closed_form) = case(&mut model);
        let (seated, work) = ogeom_bool::work::counted(|| {
            ogeom_fillet::fillet_edges(&mut model, &solid, &edges, radius, T).unwrap()
        });
        assert_eq!(work.faces_gathered == 0, trimmed, "{name}: trimmed");
        assert!(
            ogeom_algo::check(&model, &seated.shape, T)
                .unwrap()
                .is_valid(),
            "{name}: the rounding is valid"
        );
        let got = volume(&model, &seated.shape, 1e-3);
        let exact = closed_form.expect("every case has an independent volume");
        // The independent formulas and quadratures are much more accurate
        // than the 1e-3 chord; its meshing error is allotted 1e-7 relative.
        assert!(
            (got - exact).abs() <= exact * 1e-7,
            "{name}: volume {got} against the independent value {exact}"
        );
    }
}

/// A concave straight band's wedge along the edge from `corner` along `d`
/// over `length`, between the faces leaving that edge along `a` and `b`
/// (unit, square to each other and to `d`): the section between the two
/// legs of length `r` and the quarter circle of radius `r` about
/// `corner + r (a + b)` swept along the edge, its five faces sewn.
fn concave_wedge(
    model: &mut ogeom_topo::Model,
    corner: Point,
    (a, b): (Vector, Vector),
    d: Vector,
    length: f64,
    r: f64,
) -> ogeom_topo::Shape {
    use ogeom_topo::Shape;
    // `a`, `b`, `d` taken right-handed.
    let (a, b) = if a.cross(b).dot(d) > 0.0 {
        (a, b)
    } else {
        (b, a)
    };
    let dir = |v: Vector| Direction::new(v, T).unwrap();
    let ends = [corner, corner + d * length];
    // At each end: the corner, the contact on `a`'s face, on `b`'s face.
    let [near, far] = ends.map(|c| {
        [c, c + a * r, c + b * r].map(|at| (ogeom_algo::make_vertex(model, at).shape, at))
    });
    let line =
        |model: &mut ogeom_topo::Model, (from, s): &(Shape, Point), (to, e): &(Shape, Point)| {
            let curve = ogeom_geom::Curve::Line(ogeom_geom::LineCurve::segment(*s, *e, T).unwrap());
            let domain = ogeom_geom::Curve3d::domain(&curve);
            ogeom_algo::make_edge_between(model, curve, domain, from, to, T)
                .unwrap()
                .shape
        };
    // Each end's quarter circle from the contact on `a`'s face to the one
    // on `b`'s, about the axis's point there.
    let arc = |model: &mut ogeom_topo::Model, end: &[(Shape, Point); 3]| {
        let frame = Frame::new(end[0].1 + (a + b) * r, dir(-d), dir(-b), T).unwrap();
        let curve = ogeom_geom::Curve::Circle(CircleCurve::new(Circle::new(frame, r, T).unwrap()));
        let quarter = (0.0, core::f64::consts::FRAC_PI_2);
        ogeom_algo::make_edge_between(model, curve, quarter, &end[1].0, &end[2].0, T)
            .unwrap()
            .shape
    };
    let rails: Vec<Shape> = (0..3).map(|k| line(model, &near[k], &far[k])).collect();
    let [[ca, arc0, bc], [ca1, arc1, bc1]] = [&near, &far].map(|end| {
        [
            line(model, &end[0], &end[1]),
            arc(model, end),
            line(model, &end[2], &end[0]),
        ]
    });
    // Each ring counter-clockwise about its surface's normal; sewing turns
    // the faces to agree and the shell to face out.
    let face = |model: &mut ogeom_topo::Model, surface: ogeom_geom::SurfaceGeometry, ring| {
        ogeom_algo::make_face_with_pcurves(model, surface, &[ring], T)
            .unwrap()
            .shape
    };
    let plane = |at: Point, normal: Vector| {
        ogeom_geom::PlaneSurface::over(
            Plane::through(at, dir(normal)),
            (-100.0, 100.0),
            (-100.0, 100.0),
        )
        .unwrap()
        .into()
    };
    let frame = Frame::new(corner + (a + b) * r, dir(d), dir(-b), T).unwrap();
    let cylinder = ogeom_geom::SurfaceGeometry::Cylinder(
        ogeom_geom::CylinderSurface::new(
            ogeom_math::Cylinder::new(frame, r, T).unwrap(),
            (-1.0, length + 1.0),
        )
        .unwrap(),
    );
    // The cylinder's normal points away from its axis.
    let blend = vec![
        rails[1].clone(),
        arc1.clone(),
        rails[2].reversed(),
        arc0.reversed(),
    ];
    let faces = [
        face(model, plane(ends[0], d), vec![ca.clone(), arc0, bc.clone()]),
        face(
            model,
            plane(ends[1], d),
            vec![ca1.clone(), arc1, bc1.clone()],
        ),
        face(
            model,
            plane(corner, -b),
            vec![ca, rails[1].clone(), ca1.reversed(), rails[0].reversed()],
        ),
        face(
            model,
            plane(corner, -a),
            vec![rails[0].clone(), bc1.reversed(), rails[2].reversed(), bc],
        ),
        face(model, cylinder, blend),
    ];
    let sewn = ogeom_algo::sew(model, &faces, T).unwrap();
    let [shell] = sewn.shells.as_slice() else {
        panic!("the wedge sews into one shell");
    };
    let wedge = ogeom_algo::make_solid(model, std::slice::from_ref(shell))
        .unwrap()
        .shape;
    assert!(ogeom_algo::check(model, &wedge, T).unwrap().is_valid());
    wedge
}

/// A concave straight band's wedge fused on where it runs, at a pocket's
/// corner, into a concave band of its radius already filling the floor's
/// edge along the wall it meets is trimmed, cut back at the mitre where the
/// two blends cross, to the solid and history the general fuse gives. The
/// pocket, 12 x 8 and 5 deep in a 20 x 16 x 10 block, has its floor's edge
/// along one long wall, or both, filled first; the wedge then fills the
/// edge along a short wall, built from either end, on either short wall, and
/// against a first band swept either way along its edge. Each band adds
/// `(1 - pi/4) r^2` along its edge and the wedge gives back, at each corner
/// it shares, the cell past the mitre, `2 r^3 (5/6 - pi/4)`. A wedge of
/// another radius than the band it runs into is left to the general fuse.
#[test]
fn concave_bands_running_into_a_filled_corner_are_trimmed_to_the_general_fuse_s_solid() {
    use core::f64::consts::FRAC_PI_4;
    let section = 1.0 - FRAC_PI_4;
    let cell = 2.0 * (5.0 / 6.0 - FRAC_PI_4);
    let (x, y, z) = (
        Vector::new(1.0, 0.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
        Vector::new(0.0, 0.0, 1.0),
    );
    // A band along the floor's edge on the wall at `y = 4` or `y = 12`,
    // swept along `+x` or `-x`.
    let long = |m: &mut ogeom_topo::Model, wall: f64, back: bool, r: f64| {
        let leg = if wall < 8.0 { y } else { -y };
        let (from, d) = if back { (16.0, -x) } else { (4.0, x) };
        concave_wedge(m, Point::new(from, wall, 5.0), (leg, z), d, 12.0, r)
    };
    // A wedge along the floor's edge on the wall at `x = 4` or `x = 16`,
    // swept along `+y` or `-y`.
    let short = |m: &mut ogeom_topo::Model, wall: f64, back: bool, r: f64| {
        let leg = if wall < 10.0 { x } else { -x };
        let (from, d) = if back { (12.0, -y) } else { (4.0, y) };
        concave_wedge(m, Point::new(wall, from, 5.0), (leg, z), d, 8.0, r)
    };
    struct Case {
        name: &'static str,
        bands: &'static [(f64, bool)],
        wedge: (f64, bool, f64),
        trimmed: bool,
        closed_form: Option<f64>,
    }
    let one = 2720.0 + 20.0 * section - cell;
    // The unequal bands overlap at the corner over 0 <= z <= 0.5. At each
    // height their widths are the distances from the walls to their arcs;
    // the overlap's horizontal section is the product of those widths.
    let mixed_overlap = gauss(
        |theta| {
            // z = sin²(theta)/2 smooths both radicals at z=0.
            let z = 0.5 * theta.sin().powi(2);
            let long = 1.0 - (2.0 * z - z * z).sqrt();
            let short = 0.5 - (z - z * z).sqrt();
            long * short * theta.sin() * theta.cos()
        },
        0.0,
        core::f64::consts::FRAC_PI_2,
    );
    let mixed = 2720.0 + 12.0 * section + 2.0 * section - mixed_overlap;
    let cases = [
        Case {
            name: "a wedge running into the band at its far end",
            bands: &[(12.0, false)],
            wedge: (16.0, false, 1.0),
            trimmed: true,
            closed_form: Some(one),
        },
        Case {
            name: "a wedge running into the band at its near end",
            bands: &[(12.0, false)],
            wedge: (16.0, true, 1.0),
            trimmed: true,
            closed_form: Some(one),
        },
        Case {
            name: "a wedge on the other short wall",
            bands: &[(12.0, false)],
            wedge: (4.0, false, 1.0),
            trimmed: true,
            closed_form: Some(one),
        },
        Case {
            name: "a band swept the other way",
            bands: &[(12.0, true)],
            wedge: (16.0, false, 1.0),
            trimmed: true,
            closed_form: Some(one),
        },
        Case {
            name: "a wedge running into a band at both ends",
            bands: &[(4.0, false), (12.0, true)],
            wedge: (4.0, true, 1.0),
            trimmed: true,
            closed_form: Some(2720.0 + 32.0 * section - 2.0 * cell),
        },
        Case {
            name: "a wedge of a smaller radius than the band",
            bands: &[(12.0, false)],
            wedge: (16.0, false, 0.5),
            trimmed: false,
            closed_form: Some(mixed),
        },
    ];
    for case in cases {
        let name = case.name;
        let mut model = ogeom_topo::Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 16.0, 10.0), T)
            .unwrap()
            .shape;
        let at = Frame::new(Point::new(4.0, 4.0, 5.0), Direction::Z, Direction::X, T).unwrap();
        let pocket = ogeom_algo::make_box(&mut model, at, (12.0, 8.0, 6.0), T)
            .unwrap()
            .shape;
        let mut solid = ogeom_bool::cut(&mut model, &block, &pocket, T)
            .unwrap()
            .shape;
        for &(wall, back) in case.bands {
            let band = long(&mut model, wall, back, 1.0);
            solid = ogeom_bool::fuse_seated(&mut model, &solid, &band, T)
                .unwrap()
                .shape;
        }
        let (wall, back, r) = case.wedge;
        let wedge = short(&mut model, wall, back, r);
        let (seated, work) = ogeom_bool::work::counted(|| {
            ogeom_bool::fuse_seated(&mut model, &solid, &wedge, T).unwrap()
        });
        assert_eq!(
            (work.seated_fills, work.faces_gathered == 0),
            (usize::from(case.trimmed), case.trimmed),
            "{name}: filled by trimming"
        );
        assert!(
            ogeom_algo::check(&model, &seated.shape, T)
                .unwrap()
                .is_valid(),
            "{name}: the fill is valid"
        );
        let got = volume(&model, &seated.shape, 1e-3);
        let exact = case
            .closed_form
            .expect("every fill has an independent volume");
        // The section product is integrated to much better than the 1e-7
        // relative allowance for tessellation at the 1e-3 chord.
        assert!(
            (got - exact).abs() <= exact * 1e-7,
            "{name}: volume {got} against the independent value {exact}"
        );
    }
}

/// A tangent chain of straight edges and arcs takes the sum of its
/// straight spandrels and revolved circular sectors, less overlapping
/// material. A blind drill removes only its portion inside the rounded
/// profile. The cases are a block's
/// top loop once its upright edges are rounded, a stadium slot's top loop
/// (the arcs on a concave wall) and a stadium prism's top loop (its placed
/// top baked), all now built by the engine's trimmed assembly, nothing
/// trimmed; for a union a rounded pocket's floor loop; and the block's top
/// loop with a hole coming up into the sector a corner's wedge lies in, its
/// upright edges rounded each on a cylinder spanning its strip alone, now
/// the engine's non-local general route (the hole off a copy, the corner
/// built there), nothing trimmed.
#[test]
fn rims_pieces_in_a_tangent_chain_are_trimmed_to_the_general_fuse_s_solid() {
    use core::f64::consts::{FRAC_PI_4, PI};
    type Made = (ogeom_topo::Shape, Vec<ogeom_topo::Shape>, f64, Option<f64>);
    type Case<'a> = &'a dyn Fn(&mut ogeom_topo::Model) -> Made;
    let section = 1.0 - FRAC_PI_4;
    let block = |model: &mut ogeom_topo::Model, at: Point, size: (f64, f64, f64)| {
        let frame = Frame::new(at, Direction::Z, Direction::X, T).unwrap();
        ogeom_algo::make_box(model, frame, size, T).unwrap().shape
    };
    let drum = |model: &mut ogeom_topo::Model, at: Point, r: f64, h: f64| {
        let frame = Frame::new(at, Direction::Z, Direction::X, T).unwrap();
        ogeom_algo::make_cylinder(model, frame, r, h, T)
            .unwrap()
            .shape
    };
    // The edges of `shape` at height `z`, and those running upright.
    let ends = |model: &ogeom_topo::Model, e: &ogeom_topo::Shape| {
        let data = model.node(e).unwrap().data().as_edge().unwrap();
        let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            unreachable!("an edge of a solid made here has a curve");
        };
        let curve = model.geometry().curve(*curve).unwrap();
        [range.0, f64::midpoint(range.0, range.1), range.1]
            .map(|t| ogeom_geom::Curve3d::point_at(curve, t, T).unwrap())
    };
    let at_height = |model: &ogeom_topo::Model, shape: &ogeom_topo::Shape, z: f64| {
        explore_unique(model, shape, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| ends(model, e).iter().all(|p| (p.z - z).abs() < 1e-9))
            .collect::<Vec<_>>()
    };
    let upright = |model: &ogeom_topo::Model, shape: &ogeom_topo::Shape| {
        explore_unique(model, shape, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                let [p, _, q] = ends(model, e);
                (p.x - q.x).abs() < 1e-9 && (p.y - q.y).abs() < 1e-9
            })
            .collect::<Vec<_>>()
    };
    // A 20 x 20 x 10 block with its upright edges rounded at 3, and a
    // blind hole drilled up into it at `hole` where one is given, coming up
    // to 0.3 under its top; its top loop.
    let rounded_block = |model: &mut ogeom_topo::Model, hole: Option<Point>| -> Made {
        let solid = block(model, Point::ORIGIN, (20.0, 20.0, 10.0));
        let edges = upright(model, &solid);
        let mut solid = ogeom_fillet::fillet_edges(model, &solid, &edges, 3.0, T)
            .unwrap()
            .shape;
        let mut want = Some(
            4000.0
                - 4.0 * section * 9.0 * 10.0
                - 4.0 * 14.0 * section
                - rim_removed(2.0, 1.0, true),
        );
        if let Some(at) = hole {
            let tool = drum(model, at, 0.2, 9.7 - at.z);
            solid = ogeom_bool::cut(model, &solid, &tool, T).unwrap().shape;
            let plain = want.expect("the rounded block has a closed-form volume");
            // At the highest drill section, the smallest rightward reach
            // of the corner's circular wall is 17 + sqrt((2+sqrt(.51))²
            // - 1.7²) > 18.7, beyond the drill's rightmost point. At all
            // lower sections it reaches farther. The full disk is inside.
            let drilled = 9.7 * PI * 0.04;
            want = Some(plain - drilled);
        }
        let edges = at_height(model, &solid, 10.0);
        (solid, edges, 1.0, want)
    };
    let stadium_top = |model: &mut ogeom_topo::Model| -> Made {
        let solid = stadium(model, 10.0, 5.0, 4.0);
        let edges: Vec<_> = explore_unique(model, &solid, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                explore(model, e, Filter::OfType(ShapeType::Vertex))
                    .unwrap()
                    .iter()
                    .all(|v| {
                        let p = model.node(v).unwrap().data().as_vertex().unwrap().point;
                        let placed = v.transform(model.datums()).unwrap().apply(p);
                        (placed.z - 4.0).abs() < 1e-9
                    })
            })
            .collect();
        let want =
            (10.0 * 10.0 + PI * 25.0) * 4.0 - 2.0 * 10.0 * section - rim_removed(4.0, 1.0, true);
        (solid, edges, 1.0, Some(want))
    };
    // A plate with a stadium slot through it, its top loop.
    let slot = |model: &mut ogeom_topo::Model| -> Made {
        let plate = block(model, Point::new(-8.0, -6.0, -5.0), (26.0, 20.0, 10.0));
        let mut tool = block(model, Point::ORIGIN, (10.0, 8.0, 6.0));
        for x in [0.0, 10.0] {
            let end = drum(model, Point::new(x, 4.0, 0.0), 4.0, 6.0);
            tool = ogeom_bool::fuse(model, &tool, &end, T).unwrap().shape;
        }
        let solid = ogeom_bool::cut(model, &plate, &tool, T).unwrap().shape;
        let edges = at_height(model, &solid, 5.0)
            .into_iter()
            .filter(|e| {
                let [_, m, _] = ends(model, e);
                m.x > -4.5 && m.x < 14.5 && m.y > -0.5 && m.y < 8.5
            })
            .collect();
        let want = 26.0 * 20.0 * 10.0
            - (10.0 * 8.0 + PI * 16.0) * 5.0
            - 2.0 * 10.0 * section
            - rim_removed(5.0, 1.0, false);
        (solid, edges, 1.0, Some(want))
    };
    // A 10 x 10 pocket 5 deep in a 20 x 20 x 10 block, its upright edges
    // rounded at 2; its floor loop.
    let pocket = |model: &mut ogeom_topo::Model| -> Made {
        let solid = block(model, Point::ORIGIN, (20.0, 20.0, 10.0));
        let tool = block(model, Point::new(5.0, 5.0, 5.0), (10.0, 10.0, 6.0));
        let solid = ogeom_bool::cut(model, &solid, &tool, T).unwrap().shape;
        let edges: Vec<_> = upright(model, &solid)
            .into_iter()
            .filter(|e| {
                let [p, _, _] = ends(model, e);
                (4.0..16.0).contains(&p.x) && (4.0..16.0).contains(&p.y)
            })
            .collect();
        let solid = ogeom_fillet::fillet_edges(model, &solid, &edges, 2.0, T)
            .unwrap()
            .shape;
        let edges = at_height(model, &solid, 5.0);
        let want = 4000.0 - 500.0
            + 4.0 * section * 4.0 * 5.0
            + 4.0 * 6.0 * section
            + rim_removed(1.0, 1.0, true);
        (solid, edges, 1.0, Some(want))
    };
    let cases: [(&str, Case, bool); 5] = [
        (
            "a rounded block's top loop",
            &|m| rounded_block(m, None),
            false,
        ),
        ("a stadium prism's top loop", &stadium_top, false),
        ("a stadium slot's top loop", &slot, false),
        ("a rounded pocket's floor loop", &pocket, true),
        (
            "a hole into a corner's sector",
            &|m| rounded_block(m, Some(Point::new(18.5, 1.5, -1.0))),
            false,
        ),
    ];
    for (name, case, trimmed) in cases {
        let mut model = ogeom_topo::Model::new();
        let (solid, edges, radius, closed_form) = case(&mut model);
        let (seated, work) = ogeom_bool::work::counted(|| {
            ogeom_fillet::fillet_edges(&mut model, &solid, &edges, radius, T).unwrap()
        });
        assert_eq!(
            work.seated_cuts + work.seated_fills > 0,
            trimmed,
            "{name}: trimmed"
        );
        assert!(
            ogeom_algo::check(&model, &seated.shape, T)
                .unwrap()
                .is_valid(),
            "{name}: the rounding is valid"
        );
        let got = volume(&model, &seated.shape, 1e-3);
        let exact = closed_form.expect("every case has an independent volume");
        // The independent formulas and quadratures are much more accurate
        // than the 1e-3 chord; its meshing error is allotted 1e-7 relative.
        assert!(
            (got - exact).abs() <= exact * 1e-7,
            "{name}: volume {got} against the independent value {exact}"
        );
    }
}
