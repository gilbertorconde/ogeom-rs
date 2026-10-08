//! Mass properties of faces the closed-form rectangles and discs do not
//! cover: an elliptic wall, a plane bounded by an ellipse or by half of
//! one. Each is integrated round its own chart boundary, so the answer
//! does not depend on the deflection passed in.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    make_edge, make_edge_between, make_face_with_pcurves, make_prism, surface_properties,
    volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{Curve, Curve3d, EllipseCurve, LineCurve, PlaneSurface};
use ogeom::math::{Ellipse, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, VertexData};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn ellipse(a: f64, b: f64) -> Curve {
    Curve::Ellipse(EllipseCurve::new(
        Ellipse::new(Frame::WORLD, a, b, T).unwrap(),
    ))
}

fn prism_of(model: &mut Model, boundary: Vec<Shape>, height: f64) -> Shape {
    let plane = PlaneSurface::over(Plane::XY, (-20.0, 20.0), (-20.0, 20.0)).unwrap();
    let face = make_face_with_pcurves(model, plane.into(), &[boundary], T)
        .unwrap()
        .shape;
    make_prism(model, &face, Vector::new(0.0, 0.0, height), T)
        .unwrap()
        .shape
}

/// The perimeter of an ellipse, by adaptive quadrature of its speed.
fn perimeter(a: f64, b: f64) -> f64 {
    ogeom::math::integrate(
        |t: f64| (a * t.sin()).hypot(b * t.cos()),
        0.0,
        2.0 * PI,
        1e-13,
    )
    .unwrap()
}

/// A 20 by 10 ellipse padded 4 high measures to rounding at the default
/// deflection, and says it was not meshed.
#[test]
fn an_elliptic_prism_measures_within_a_part_in_a_thousand() {
    let mut model = Model::new();
    let curve = ellipse(10.0, 5.0);
    let range = curve.domain();
    let edge = make_edge(&mut model, curve, range, T).unwrap().shape;
    let prism = prism_of(&mut model, vec![edge], 4.0);

    let v = volume_properties(&model, &prism, Deflection::default(), T).unwrap();
    let exact = PI * 10.0 * 5.0 * 4.0;
    assert!((v.mass - exact).abs() < 1e-9 * exact, "{}", v.mass);
    assert_eq!(v.deflection, 0.0, "integrated, not meshed");
    let centre = v.centre;
    assert!(
        centre.distance(Point::new(0.0, 0.0, 2.0)) < 1e-9,
        "{centre:?}"
    );

    let s = surface_properties(&model, &prism, Deflection::default(), T).unwrap();
    let area = 2.0 * PI * 10.0 * 5.0 + perimeter(10.0, 5.0) * 4.0;
    assert!(
        (s.mass - area).abs() < 1e-9 * area,
        "{} against {area}",
        s.mass
    );
}

/// Half of that ellipse, closed along its major axis: an arc and a line
/// on the plane, an elliptic wall and a flat one round the side.
#[test]
fn half_an_elliptic_prism_measures_half() {
    let mut model = Model::new();
    let (east, west) = (Point::new(10.0, 0.0, 0.0), Point::new(-10.0, 0.0, 0.0));
    let (e, w) = (
        model.add_vertex(VertexData::new(east)),
        model.add_vertex(VertexData::new(west)),
    );
    let arc = make_edge_between(&mut model, ellipse(10.0, 5.0), (0.0, PI), &e, &w, T)
        .unwrap()
        .shape;
    let chord = LineCurve::segment(west, east, T).unwrap();
    let range = chord.domain();
    let base = make_edge_between(&mut model, Curve::Line(chord), range, &w, &e, T)
        .unwrap()
        .shape;
    let half = prism_of(&mut model, vec![arc, base], 4.0);

    let v = volume_properties(&model, &half, Deflection::default(), T).unwrap();
    let exact = PI * 10.0 * 5.0 * 4.0 / 2.0;
    assert!((v.mass - exact).abs() < 1e-9 * exact, "{}", v.mass);
    assert_eq!(v.deflection, 0.0);
    // The centroid of a half ellipse stands 4b/(3 pi) off its major axis.
    let centre = v.centre;
    let y = 4.0 * 5.0 / (3.0 * PI);
    assert!(
        centre.distance(Point::new(0.0, y, 2.0)) < 1e-9,
        "{centre:?}"
    );
}

/// Half the elliptic prism and the three-quarter drum below, each placed
/// by similarities that scale, turn, mirror and move: their faces are
/// still integrated round their chart boundaries, and volume, area and
/// centre follow the placement to rounding.
#[test]
fn chart_faces_under_a_scaling_placement_measure_exactly() {
    use ogeom::math::{Axis, Direction, Transform};
    let mut model = Model::new();
    let (east, west) = (Point::new(10.0, 0.0, 0.0), Point::new(-10.0, 0.0, 0.0));
    let (e, w) = (
        model.add_vertex(VertexData::new(east)),
        model.add_vertex(VertexData::new(west)),
    );
    let arc = make_edge_between(&mut model, ellipse(10.0, 5.0), (0.0, PI), &e, &w, T)
        .unwrap()
        .shape;
    let chord = LineCurve::segment(west, east, T).unwrap();
    let range = chord.domain();
    let base = make_edge_between(&mut model, Curve::Line(chord), range, &w, &e, T)
        .unwrap()
        .shape;
    let half = prism_of(&mut model, vec![arc, base], 4.0);

    let below = Frame::new(Point::new(0.0, 0.0, -3.0), Direction::Z, Direction::X, T).unwrap();
    let drum = ogeom::algo::make_cylinder(&mut model, below, 5.0, 6.0, T)
        .unwrap()
        .shape;
    let corner = Frame::new(Point::new(0.0, 0.0, -20.0), Direction::Z, Direction::X, T).unwrap();
    let quarter = ogeom::algo::make_box(&mut model, corner, (20.0, 20.0, 40.0), T)
        .unwrap()
        .shape;
    let rest = ogeom::boolean::cut(&mut model, &drum, &quarter, T)
        .unwrap()
        .shape;

    // Volume, area and centre of each as it stands unplaced.
    let lean = 4.0 * 5.0 / (3.0 * PI);
    let off = -4.0 * 5.0 / (9.0 * PI);
    let shapes = [
        (
            &half,
            PI * 10.0 * 5.0 * 4.0 / 2.0,
            PI * 10.0 * 5.0 + perimeter(10.0, 5.0) * 2.0 + 20.0 * 4.0,
            Point::new(0.0, lean, 2.0),
        ),
        (
            &rest,
            0.75 * PI * 25.0 * 6.0,
            1.5 * PI * 25.0 + 0.75 * 2.0 * PI * 5.0 * 6.0 + 2.0 * 5.0 * 6.0,
            Point::new(off, off, 0.0),
        ),
    ];
    let placements = [
        Transform::scaling(Point::ORIGIN, 0.5, T).unwrap(),
        Transform::translation(Vector::new(7.0, -3.0, 11.0))
            * Transform::rotation(Axis::new(Point::ORIGIN, Direction::X), 0.7)
            * Transform::scaling(Point::new(1.0, 2.0, 3.0), 2.0, T).unwrap(),
        Transform::plane_mirror(Point::new(0.0, 4.0, 0.0), Direction::Y)
            * Transform::scaling(Point::new(-2.0, 0.0, 1.0), 0.25, T).unwrap(),
        Transform::scaling(Point::new(2.0, 0.0, 0.0), -3.0, T).unwrap(),
    ];
    for (shape, volume, area, centre) in shapes {
        for placement in placements {
            let placed = ogeom::algo::transformed(&mut model, shape, placement)
                .unwrap()
                .shape;
            let s = placement.scale_factor().abs();
            let at = placement.apply(centre);
            let v = volume_properties(&model, &placed, Deflection::default(), T).unwrap();
            let a = surface_properties(&model, &placed, Deflection::default(), T).unwrap();
            for (found, exact) in [(v, volume * s.powi(3)), (a, area * s * s)] {
                assert_eq!(found.deflection, 0.0, "integrated, not meshed");
                assert!(
                    (found.mass - exact).abs() < 1e-9 * exact,
                    "{} against {exact}",
                    found.mass
                );
            }
            assert!(
                v.centre.distance(at) < 1e-9 * (1.0 + at.to_vector().magnitude()),
                "{:?} against {at:?}",
                v.centre
            );
        }
    }
}

/// A drum with a quarter cut away, along its seam and across its axis:
/// its caps are concave, three quarters of a disc, and its wall keeps the
/// seam down one column only. It measures exactly, which needs every face
/// read the right way out.
#[test]
fn a_three_quarter_drum_measures_exactly() {
    use ogeom::math::Direction;
    let mut model = Model::new();
    let below = Frame::new(Point::new(0.0, 0.0, -3.0), Direction::Z, Direction::X, T).unwrap();
    let drum = ogeom::algo::make_cylinder(&mut model, below, 5.0, 6.0, T)
        .unwrap()
        .shape;
    let corner = Frame::new(Point::new(0.0, 0.0, -20.0), Direction::Z, Direction::X, T).unwrap();
    let quarter = ogeom::algo::make_box(&mut model, corner, (20.0, 20.0, 40.0), T)
        .unwrap()
        .shape;
    let rest = ogeom::boolean::cut(&mut model, &drum, &quarter, T)
        .unwrap()
        .shape;
    let v = volume_properties(&model, &rest, Deflection::default(), T).unwrap();
    let exact = 0.75 * PI * 25.0 * 6.0;
    assert_eq!(v.deflection, 0.0, "integrated, not meshed");
    assert!((v.mass - exact).abs() < 1e-9 * exact, "{}", v.mass);
}

/// The area of the cone wall of a frustum of radius 6 at its base and 3 at
/// its top, 10 high, that a drill of radius `r` lying along y on its base
/// through its rim takes away: the part of the wall within `r` of the line
/// `x = 6, z = r`. At height `z` the wall's circle of radius `rho` is inside
/// the drill where `rho cos(theta)` exceeds `6 - s`, `s` the drill's half
/// width there, up to the height `2r / 1.09` where the drill leaves the
/// wall. Its square-root ends are smoothed by a substitution at each.
fn drilled_cone_wall(r: f64) -> f64 {
    let slant = 1.09_f64.sqrt();
    let width = |z: f64| {
        let rho = 6.0 - 0.3 * z;
        let s = (2.0 * r * z - z * z).max(0.0).sqrt();
        2.0 * ((6.0 - s) / rho).clamp(-1.0, 1.0).acos() * rho * slant
    };
    let top = 2.0 * r / 1.09;
    let half = 0.5 * top;
    let low = ogeom::math::integrate(
        |t| width(half * t.powi(4)) * 4.0 * half * t.powi(3),
        0.0,
        1.0,
        1e-13,
    )
    .unwrap();
    let high = ogeom::math::integrate(
        |t| width(top - half * t * t) * 2.0 * half * t,
        0.0,
        1.0,
        1e-13,
    )
    .unwrap();
    low + high
}

/// A frustum drilled by a drill lying on its base, the drill's wall
/// touching the base's circle at the cone wall's seam. The section on the
/// wall is fitted, its pieces' ends a few hundred-thousandths off the
/// vertices they meet, and the cut's wall wraps the whole chart round the
/// hole. The cut's wall measures the lateral area less the hole's, the
/// common's wall the hole's, and cut and common add up to the frustum to
/// a part in a billion.
#[test]
fn a_frustum_drilled_through_its_rim_measures_its_hole() {
    use ogeom::algo::{make_cone, make_cylinder};
    use ogeom::math::Direction;
    use ogeom::topo::{NodeData, ShapeType, explore_unique};
    let fine = Deflection::with_chord(1e-3).unwrap();
    let cone_walls = |model: &Model, shape: &Shape| -> f64 {
        explore_unique(model, shape, ShapeType::Face)
            .unwrap()
            .iter()
            .filter(|face| {
                let NodeData::Face(data) = model.node(face).unwrap().data() else {
                    return false;
                };
                matches!(
                    model.geometry().surface(data.surface),
                    Some(ogeom::geom::SurfaceGeometry::Cone(_))
                )
            })
            .map(|face| surface_properties(model, face, fine, T).unwrap().mass)
            .sum()
    };
    let lateral = PI * 9.0 * 109.0_f64.sqrt();
    let frustum = PI * 10.0 / 3.0 * (36.0 + 18.0 + 9.0);
    for r in [0.3, 0.5] {
        let mut model = Model::new();
        let part = make_cone(&mut model, Frame::WORLD, 6.0, 3.0, 10.0, T)
            .unwrap()
            .shape;
        let frame = Frame::new(
            Point::new(6.0, -8.0, r),
            Direction::Y,
            Direction::new(Vector::new(0.0, 0.0, -1.0), T).unwrap(),
            T,
        )
        .unwrap();
        let drill = make_cylinder(&mut model, frame, r, 16.0, T).unwrap().shape;
        let cut = ogeom::boolean::cut(&mut model, &part, &drill, T)
            .unwrap()
            .shape;
        let common = ogeom::boolean::common(&mut model, &part, &drill, T)
            .unwrap()
            .shape;

        let hole = drilled_cone_wall(r);
        let (kept, taken) = (cone_walls(&model, &cut), cone_walls(&model, &common));
        // The two walls share the section's pcurves, so between them they
        // cover the lateral area to rounding.
        assert!(
            (kept + taken - lateral).abs() < 1e-9 * lateral,
            "r {r}: the walls {kept} + {taken} against {lateral}"
        );
        // Each against the hole: the section is fitted to within 2e-5 of
        // the true one along its 2 mm or so, which moves the hole's area by
        // a few millionths.
        assert!(
            (kept - (lateral - hole)).abs() < 5e-6,
            "r {r}: the cut's wall {kept} against {}",
            lateral - hole
        );
        assert!(
            (taken - hole).abs() < 5e-6,
            "r {r}: the common's wall {taken} against {hole}"
        );

        let (v_cut, v_common) = (
            volume_properties(&model, &cut, fine, T).unwrap(),
            volume_properties(&model, &common, fine, T).unwrap(),
        );
        assert_eq!(v_cut.deflection, 0.0, "integrated, not meshed");
        assert_eq!(v_common.deflection, 0.0, "integrated, not meshed");
        let sum = v_cut.mass + v_common.mass;
        assert!(
            (sum - frustum).abs() < 1e-9 * frustum,
            "r {r}: cut {} + common {} against {frustum}",
            v_cut.mass,
            v_common.mass
        );
    }
}

/// The common of a torus of radii 10 and 3 and a drill along its axis of
/// radius `r` centred at `(x, y)`: the torus's thickness along z, twice the
/// root of 9 less the square of the distance from the tube's centre
/// circle, integrated over the drill's disc. Taken ring by ring round the
/// torus's axis, each ring's angle inside the disc in closed form; the
/// square-root ends at the tube's walls and where the rings start to leave
/// the disc are smoothed by a sine substitution on each stretch.
fn torus_core(x: f64, y: f64, r: f64) -> f64 {
    let d = x.hypot(y);
    let f = |rho: f64| {
        let left = 9.0 - (rho - 10.0) * (rho - 10.0);
        let angle = if rho + d <= r {
            2.0 * PI
        } else if rho >= d + r || rho <= d - r {
            0.0
        } else {
            2.0 * ((rho * rho + d * d - r * r) / (2.0 * rho * d))
                .clamp(-1.0, 1.0)
                .acos()
        };
        2.0 * left.max(0.0).sqrt() * rho * angle
    };
    let mut breaks = vec![7.0, 13.0];
    breaks.extend(
        [d - r, d + r, r - d]
            .into_iter()
            .filter(|b| *b > 7.0 && *b < 13.0),
    );
    breaks.sort_by(f64::total_cmp);
    breaks
        .windows(2)
        .map(|pair| {
            let (mid, half) = (0.5 * (pair[0] + pair[1]), 0.5 * (pair[1] - pair[0]));
            ogeom::math::integrate(
                |phi| f(mid + half * phi.sin()) * half * phi.cos(),
                -0.5 * PI,
                0.5 * PI,
                1e-13,
            )
            .unwrap()
        })
        .sum()
}

/// A torus drilled along its axis by drills whose wall passes through the
/// vertex where its seams cross on the outer equator, square to the
/// equator there. Near the vertex the section hugs the seam's meridian
/// circle, and the boolean bounds the drill's wall there by the meridian:
/// the wall's pcurve, fitted to the true section, stands up to 2e-5 off
/// the edge's circle, and the fitted sections' own pcurves stand off their
/// edges by up to their stated 1e-5. The strips between them are integrated
/// with the faces, and the common measures the torus's core within the
/// drill to a part in a trillion, where leaving out the strips beside the
/// circle misses by 2.6e-8 and those beside the fitted sections by 1.5e-9.
#[test]
fn a_torus_drilled_through_its_seams_vertex_measures_its_core() {
    use ogeom::algo::{make_cylinder, make_torus};
    use ogeom::math::Direction;
    let fine = Deflection::with_chord(1e-3).unwrap();
    for (r, seam) in [
        (3.721_499_296_241_649_5, Vector::new(0.0, 1.0, 0.0)),
        (-4.243_598_220_801_017, Vector::new(-1.0, 0.0, 0.0)),
    ] {
        let mut model = Model::new();
        let part = make_torus(&mut model, Frame::WORLD, 10.0, 3.0, T)
            .unwrap()
            .shape;
        let frame = Frame::new(
            Point::new(13.0, r, -6.725_587_333_175_463),
            Direction::Z,
            Direction::new(seam, T).unwrap(),
            T,
        )
        .unwrap();
        let drill = make_cylinder(&mut model, frame, r.abs(), 13.451_174_666_350_926, T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &part, &drill, T)
            .unwrap()
            .shape;
        let common = ogeom::boolean::common(&mut model, &part, &drill, T)
            .unwrap()
            .shape;
        let (v_cut, v_common) = (
            volume_properties(&model, &cut, fine, T).unwrap(),
            volume_properties(&model, &common, fine, T).unwrap(),
        );
        assert_eq!(v_cut.deflection, 0.0, "r {r}: integrated, not meshed");
        assert_eq!(v_common.deflection, 0.0, "r {r}: integrated, not meshed");
        let core = torus_core(13.0, r, r.abs());
        assert!(
            (v_common.mass - core).abs() < 1e-12 * core,
            "r {r}: common {} against {core} integrated",
            v_common.mass
        );
        let torus = 2.0 * PI * PI * 10.0 * 9.0;
        assert!(
            (v_cut.mass + v_common.mass - torus).abs() < 1e-12 * torus,
            "r {r}: cut {} + common {} against {torus}",
            v_cut.mass,
            v_common.mass
        );
    }
}

/// A box rounded on every edge, drilled along y by a drill whose wall runs
/// along the box's side face where a round leaves it. The drill's section
/// on a corner's sphere ends at a vertex whose ball holds both pieces'
/// ends, 6e-5 apart, wider than the vertex's tolerance but within its
/// diameter. Cut and common are integrated on their surfaces, and add up
/// to the rounded box to a part in a billion.
#[test]
fn a_rounded_box_drilled_along_a_round_s_edge_measures_exactly() {
    use ogeom::algo::{make_box, make_cylinder};
    use ogeom::math::Direction;
    use ogeom::topo::{ShapeType, explore_unique};
    let fine = Deflection::with_chord(1e-3).unwrap();
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (20.0, 14.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let part = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
        .unwrap()
        .shape;
    let r = 0.700_445_137_892_257_4;
    let frame = Frame::new(
        Point::new(0.0, -12.638_181_458_367_212, 8.0 + r),
        Direction::Y,
        Direction::new(Vector::new(-1.0, 0.0, 0.0), T).unwrap(),
        T,
    )
    .unwrap();
    let drill = make_cylinder(&mut model, frame, r, 29.276_362_916_734_424, T)
        .unwrap()
        .shape;
    let whole = volume_properties(&model, &part, fine, T).unwrap();
    let cut = ogeom::boolean::cut(&mut model, &part, &drill, T)
        .unwrap()
        .shape;
    let common = ogeom::boolean::common(&mut model, &part, &drill, T)
        .unwrap()
        .shape;
    let (v_cut, v_common) = (
        volume_properties(&model, &cut, fine, T).unwrap(),
        volume_properties(&model, &common, fine, T).unwrap(),
    );
    for (name, v) in [("part", &whole), ("cut", &v_cut), ("common", &v_common)] {
        assert_eq!(v.deflection, 0.0, "the {name} is integrated, not meshed");
    }
    assert!(
        (v_cut.mass + v_common.mass - whole.mass).abs() < 1e-9 * whole.mass,
        "cut {} + common {} against {}",
        v_cut.mass,
        v_common.mass,
        whole.mass
    );
}

/// A prism over a sector of radius `r` and angle `angle`, `h` high, whose
/// arc edge (one edge, at the bottom and moved to the top) takes the arc's
/// straight chord for its curve and states the chord's sag as its
/// tolerance. The round wall, a chart rectangle, keeps the arc for its
/// pcurves; the caps are bounded by the chord where `chord_caps`, and keep
/// the arc for theirs otherwise.
fn sector_prism_with_a_chord(r: f64, angle: f64, h: f64, chord_caps: bool) -> (Model, Shape) {
    use ogeom::core::Tolerance;
    use ogeom::geom::{CircleCurve, Curve2d as _, Line2d, PlanarCurve, SurfaceGeometry};
    use ogeom::math::Circle;
    use ogeom::topo::{EdgeRepr, NodeData, ShapeType, explore_unique};
    let mut model = Model::new();
    let (o, a, b) = (
        Point::ORIGIN,
        Point::new(r, 0.0, 0.0),
        Point::new(r * angle.cos(), r * angle.sin(), 0.0),
    );
    let [vo, va, vb] = [o, a, b].map(|p| model.add_vertex(VertexData::new(p)));
    let line = |model: &mut Model, from: &Shape, p: Point, to: &Shape, q: Point| {
        let curve = LineCurve::segment(p, q, T).unwrap();
        let range = curve.domain();
        make_edge_between(model, Curve::Line(curve), range, from, to, T)
            .unwrap()
            .shape
    };
    let out = line(&mut model, &vo, o, &va, a);
    let back = line(&mut model, &vb, b, &vo, o);
    let arc = Curve::Circle(CircleCurve::new(Circle::new(Frame::WORLD, r, T).unwrap()));
    let arc_edge = make_edge_between(&mut model, arc, (0.0, angle), &va, &vb, T)
        .unwrap()
        .shape;
    let prism = prism_of(&mut model, vec![out, arc_edge, back], h);

    let rim = explore_unique(&model, &prism, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let data = model.node(e).unwrap().data().as_edge().unwrap();
            let Some(EdgeRepr::Curve3d { curve, .. }) = data.curve3d() else {
                return false;
            };
            matches!(model.geometry().curve(*curve), Some(Curve::Circle(_)))
        })
        .expect("the prism has an arc");
    let data = model.node(&rim).unwrap().data().as_edge().unwrap().clone();
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        unreachable!()
    };
    let (start, end) = {
        let arc = model.geometry().curve(*curve).unwrap();
        (
            arc.point_at(range.0, T).unwrap(),
            arc.point_at(range.1, T).unwrap(),
        )
    };
    let chord = LineCurve::segment(start, end, T).unwrap();
    let chord_range = chord.domain();
    let chord = model.geometry_mut().add_curve(Curve::Line(chord));
    let mut representations = data.representations.clone();
    for repr in &mut representations {
        match repr {
            EdgeRepr::Curve3d { curve, range, .. } => {
                *curve = chord;
                *range = chord_range;
            }
            EdgeRepr::PCurve {
                curve,
                range,
                surface,
                ..
            } if chord_caps
                && matches!(
                    model.geometry().surface(*surface),
                    Some(SurfaceGeometry::Plane(_))
                ) =>
            {
                let old = model.geometry().pcurve(*curve).unwrap().clone();
                let line = Line2d::segment(
                    old.point_at(range.0, T).unwrap(),
                    old.point_at(range.1, T).unwrap(),
                    T,
                )
                .unwrap();
                *range = line.domain();
                *curve = model.geometry_mut().add_pcurve(PlanarCurve::Line(line));
            }
            _ => {}
        }
    }
    let NodeData::Edge(edge) = model.node_mut(&rim).unwrap().data_mut() else {
        unreachable!()
    };
    edge.representations = representations.into_iter().collect();
    edge.assert_same_parameter(false);
    let sag = r * (1.0 - (angle * 0.5).cos());
    model
        .widen(&rim, Tolerance::new(sag * 1.001).unwrap())
        .unwrap();
    (model, prism)
}

/// A sector prism of radius 10 and angle 0.6, 5 high, whose arc edge takes
/// the arc's chord for its curve and states the chord's sag, 0.45 (see
/// [`sector_prism_with_a_chord`]). The sliver between chord and arc at each
/// end is closed by the strips integrated with the faces whose pcurves keep
/// the arc. With caps bounded by the chord, the wall's strip closes it:
/// the wall's area is its own plus both slivers'. With caps that keep the
/// arc too, the caps' strips lie back over them and take the slivers away,
/// and the wall's strips add them, so the two cancel in the volume. Either
/// way the volume and the whole surface are the sector prism's, and each
/// cap is the triangle under the chord, each to a part in a billion
/// against closed forms. Without the strips the volume misses the flux
/// through the slivers.
#[test]
fn a_chord_standing_off_its_arc_by_its_sag_closes_with_the_strips() {
    use ogeom::geom::SurfaceGeometry;
    use ogeom::topo::{ShapeType, explore_unique};
    let (r, angle, h) = (10.0_f64, 0.6_f64, 5.0_f64);
    let sector = 0.5 * r * r * angle;
    let sliver = 0.5 * r * r * (angle - angle.sin());
    let triangle = 0.5 * r * r * angle.sin();
    for chord_caps in [true, false] {
        let (model, prism) = sector_prism_with_a_chord(r, angle, h, chord_caps);
        let v = volume_properties(&model, &prism, Deflection::default(), T).unwrap();
        assert_eq!(v.deflection, 0.0, "{chord_caps}: integrated, not meshed");
        assert!(
            (v.mass - sector * h).abs() < 1e-9 * sector * h,
            "{chord_caps}: volume {} against {}",
            v.mass,
            sector * h
        );
        let whole = 2.0 * sector + r * angle * h + 2.0 * r * h;
        let s = surface_properties(&model, &prism, Deflection::default(), T).unwrap();
        assert_eq!(s.deflection, 0.0, "{chord_caps}: integrated, not meshed");
        assert!(
            (s.mass - whole).abs() < 1e-9 * whole,
            "{chord_caps}: area {} against {whole}",
            s.mass
        );
        for face in explore_unique(&model, &prism, ShapeType::Face).unwrap() {
            let data = model.node(&face).unwrap().data().as_face().unwrap();
            let want = match model.geometry().surface(data.surface) {
                Some(SurfaceGeometry::Plane(p)) if p.plane().normal().vector().z.abs() > 0.5 => {
                    triangle
                }
                Some(SurfaceGeometry::Plane(_)) => r * h,
                _ => r * angle * h + 2.0 * sliver,
            };
            let a = surface_properties(&model, &face, Deflection::default(), T).unwrap();
            assert!(
                (a.mass - want).abs() < 1e-9 * want,
                "{chord_caps}: face {} against {want}",
                a.mass
            );
        }
    }
}

/// A cube of side 10 whose top face's pcurve along one edge runs on the
/// edge's line from a twelfth of the way to eleven twelfths, and bends off
/// it towards each end to stand 1e-4 aside at the vertices, which state
/// that much. Every inner point the strip test samples lies on the line;
/// the ends do not, and the slit between pcurve and line is closed by the
/// strip at whichever face the moments are taken from: the volume is the
/// cube's to rounding with each face listed first, where leaving the slit
/// open moves it by a few parts in ten million.
#[test]
fn a_pcurve_straying_off_its_line_only_towards_its_ends_closes() {
    use ogeom::algo::{make_box, make_shell, make_solid};
    use ogeom::core::Tolerance;
    use ogeom::geom::{BSpline2d, Curve2d as _, PlanarCurve, SurfaceGeometry};
    use ogeom::math::KnotVector;
    use ogeom::topo::{EdgeRepr, NodeData, ShapeType, explore_unique};
    let (side, aside) = (10.0, 1e-4);
    let mut model = Model::new();
    let cube = make_box(&mut model, Frame::WORLD, (side, side, side), T)
        .unwrap()
        .shape;
    let faces = explore_unique(&model, &cube, ShapeType::Face).unwrap();
    let top = faces
        .iter()
        .find(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::Plane(p)) if p.plane().origin().z > side / 2.0
            )
        })
        .expect("the cube has a top")
        .clone();
    let surface = model.node(&top).unwrap().data().as_face().unwrap().surface;
    let edge = explore_unique(&model, &top, ShapeType::Edge).unwrap()[0].clone();
    let mut data = model.node(&edge).unwrap().data().as_edge().unwrap().clone();
    for repr in &mut data.representations {
        if let EdgeRepr::PCurve {
            curve,
            range,
            surface: on,
            ..
        } = repr
            && *on == surface
        {
            let old = model.geometry().pcurve(*curve).unwrap().clone();
            let (a, b) = (
                old.point_at(range.0, T).unwrap(),
                old.point_at(range.1, T).unwrap(),
            );
            let along = b - a;
            let off = along.perpendicular() * (aside / along.magnitude());
            let (r0, r1) = *range;
            let knots = KnotVector::new(
                vec![
                    r0,
                    r0,
                    r0 + (r1 - r0) / 12.0,
                    r0 + (r1 - r0) * 11.0 / 12.0,
                    r1,
                    r1,
                ],
                1,
            )
            .unwrap();
            let control = vec![
                a + off,
                a + along / 12.0,
                a + along * (11.0 / 12.0),
                b + off,
            ];
            let bent = BSpline2d::new(knots, control, T).unwrap();
            *curve = model.geometry_mut().add_pcurve(PlanarCurve::BSpline(bent));
        }
    }
    let NodeData::Edge(slot) = model.node_mut(&edge).unwrap().data_mut() else {
        unreachable!()
    };
    **slot = data;
    model
        .widen(&edge, Tolerance::new(aside * 1.001).unwrap())
        .unwrap();
    let cube_volume = side * side * side;
    for first in 0..faces.len() {
        let mut listed = faces.clone();
        listed.rotate_left(first);
        let shell = make_shell(&mut model, &listed).unwrap().shape;
        let solid = make_solid(&mut model, &[shell]).unwrap().shape;
        let v = volume_properties(&model, &solid, Deflection::default(), T).unwrap();
        assert_eq!(v.deflection, 0.0, "first {first}: integrated, not meshed");
        assert!(
            (v.mass - cube_volume).abs() < 1e-12 * cube_volume,
            "first {first}: volume {} against {cube_volume}",
            v.mass
        );
    }
}
