//! A many-faced closed solid assembled from nothing but vertices, curves and
//! the builders — `make_wire`, `make_face_with_pcurves`, `make_shell`,
//! `make_solid` — and then measured against a closed form.
//!
//! This is the shape of use the builders get from a caller that is *not* one
//! of the kernel's own modelling operations: dozens of faces handed over at
//! once, every wire assembled by hand from segments the caller chose, and the
//! result expected to come back a valid closed solid. A regular prism is the
//! smallest thing with that shape which also has an exact volume, so the
//! measurement is an oracle rather than a bound: a regular `n`-gon of
//! circumradius `r` has area `n·r²·sin(2π/n)/2`, and the prism is that times
//! its height, to the last bit the arithmetic allows.
//!
//! The face builder here is the one that attaches pcurves. `make_face_on`
//! takes wires and a surface and nothing else, which leaves every edge on the
//! face without a curve in the face's own parameters — `check` reports that as
//! broken and the face does not triangulate. A caller assembling a face by
//! hand wants the pcurve-attaching one; the bare builder is for callers that
//! already computed the pcurves themselves.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::{LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Plane, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A segment edge between two existing vertices.
fn segment(model: &mut Model, from: &Shape, to: &Shape, a: Point, b: Point) -> Shape {
    ogeom::algo::make_edge_between(
        model,
        LineCurve::segment(a, b, T).unwrap().into(),
        (0.0, a.distance(b)),
        from,
        to,
        T,
    )
    .unwrap()
    .shape
}

/// A planar face through `origin` with normal `normal`, bounded by `edges`.
/// `reference` fixes the chart's u direction and must not be parallel to the
/// normal — a wall's outward normal is horizontal, so `Z` serves; a cap's is
/// `Z`, so `X` does.
fn face(
    model: &mut Model,
    origin: Point,
    normal: Direction,
    reference: Direction,
    edges: Vec<Shape>,
) -> Shape {
    let frame = Frame::new(origin, normal, reference, T).unwrap();
    ogeom::algo::make_face_with_pcurves(
        model,
        SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(frame))),
        &[edges],
        T,
    )
    .unwrap()
    .shape
}

#[test]
fn a_hand_built_prism_closes_and_measures_its_closed_form() {
    const SIDES: usize = 17; // odd, so no wall is parallel to its opposite
    const RADIUS: f64 = 5.0;
    const HEIGHT: f64 = 12.0;

    let mut model = Model::new();

    // The two rings of corners, and one vertex apiece — shared between the
    // cap that uses them and the two walls that meet there, which is the
    // whole point of handing the builders vertices rather than points.
    let ring = |z: f64| -> Vec<Point> {
        (0..SIDES)
            .map(|i| {
                #[allow(clippy::cast_precision_loss, reason = "a corner index")]
                let a = core::f64::consts::TAU * (i as f64) / (SIDES as f64);
                Point::new(RADIUS * a.cos(), RADIUS * a.sin(), z)
            })
            .collect()
    };
    let (bottom_pts, top_pts) = (ring(0.0), ring(HEIGHT));
    let mk = |model: &mut Model, pts: &[Point]| -> Vec<Shape> {
        pts.iter()
            .map(|p| ogeom::algo::make_vertex(model, *p).shape)
            .collect()
    };
    let bottom_vs = mk(&mut model, &bottom_pts);
    let top_vs = mk(&mut model, &top_pts);

    // The rim edges, then the uprights. Every one is used by exactly two
    // faces, which is what makes the shell closable.
    let rim = |model: &mut Model, vs: &[Shape], pts: &[Point]| -> Vec<Shape> {
        (0..SIDES)
            .map(|i| {
                let j = (i + 1) % SIDES;
                segment(model, &vs[i], &vs[j], pts[i], pts[j])
            })
            .collect()
    };
    let bottom_edges = rim(&mut model, &bottom_vs, &bottom_pts);
    let top_edges = rim(&mut model, &top_vs, &top_pts);
    let uprights: Vec<Shape> = (0..SIDES)
        .map(|i| {
            segment(
                &mut model,
                &bottom_vs[i],
                &top_vs[i],
                bottom_pts[i],
                top_pts[i],
            )
        })
        .collect();

    let mut faces = Vec::with_capacity(SIDES + 2);
    faces.push(face(
        &mut model,
        bottom_pts[0],
        Direction::Z,
        Direction::X,
        bottom_edges.clone(),
    ));
    faces.push(face(
        &mut model,
        top_pts[0],
        Direction::Z,
        Direction::X,
        top_edges.clone(),
    ));
    for i in 0..SIDES {
        let j = (i + 1) % SIDES;
        let mid = Point::new(
            f64::midpoint(bottom_pts[i].x, bottom_pts[j].x),
            f64::midpoint(bottom_pts[i].y, bottom_pts[j].y),
            0.0,
        );
        let outward = Direction::from_coords(mid.x, mid.y, 0.0, T).unwrap();
        faces.push(face(
            &mut model,
            bottom_pts[i],
            outward,
            Direction::Z,
            // Head to tail all the way round: along the bottom rim, up the
            // far upright, back along the top rim against its own sense, and
            // down the near upright against its own. `make_wire` insists on
            // this and says which pair failed to meet — an orientation the
            // caller got wrong is a gap in every face built on the wire.
            vec![
                bottom_edges[i].clone(),
                uprights[j].clone(),
                top_edges[i].reversed(),
                uprights[i].reversed(),
            ],
        ));
    }

    let shell = ogeom::algo::make_shell(&mut model, &faces).unwrap().shape;
    assert!(
        ogeom::algo::is_shell_closed(&model, &shell).unwrap(),
        "every edge used twice, so the shell closes"
    );
    let solid = ogeom::algo::make_solid(&mut model, std::slice::from_ref(&shell))
        .unwrap()
        .shape;

    let check = ogeom::algo::check(&model, &solid, T).unwrap();
    assert!(check.is_valid(), "{check}");
    assert_eq!(
        explore_unique(&model, &solid, ShapeType::Face)
            .unwrap()
            .len(),
        SIDES + 2
    );
    assert_eq!(
        explore_unique(&model, &solid, ShapeType::Edge)
            .unwrap()
            .len(),
        3 * SIDES
    );
    assert_eq!(
        explore_unique(&model, &solid, ShapeType::Vertex)
            .unwrap()
            .len(),
        2 * SIDES
    );

    // The oracle. A regular n-gon of circumradius r has area n·r²·sin(2π/n)/2.
    #[allow(clippy::cast_precision_loss, reason = "a side count")]
    let n = SIDES as f64;
    let exact = 0.5 * n * RADIUS * RADIUS * (core::f64::consts::TAU / n).sin() * HEIGHT;
    let measured = ogeom::algo::volume_properties(&model, &solid, Deflection::default(), T)
        .unwrap()
        .mass;
    // A prism is flat everywhere, so the tessellation is exact and the
    // deflection buys nothing — the tolerance here is arithmetic, not
    // sampling.
    assert!(
        (measured - exact).abs() < 1e-9,
        "the prism measures its closed form: {measured} against {exact}"
    );
}
