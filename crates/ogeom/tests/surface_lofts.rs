//! Sheets: ruled surfaces, lofts through open and closed sections, and
//! sweeps of open profiles, measured against their closed forms.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    Spacing, check, interpolate, make_edge, make_polygon, make_wire, project_on_surface,
};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve, Curve2d as _, Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom::math::{Circle, Direction, Frame, Point, Vector};
use ogeom::offset::{
    PipeLaw, make_loft_surface, make_ruled, make_sweep_surface, make_sweep_two_rails,
};
use ogeom::topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// An arc of `radius` about `centre` in the plane square to `normal`, from
/// angle `a` to `b` measured from `x`.
fn arc(
    model: &mut Model,
    centre: Point,
    normal: Vector,
    x: Vector,
    radius: f64,
    a: f64,
    b: f64,
) -> Shape {
    let frame = Frame::new(
        centre,
        Direction::new(normal, T).unwrap(),
        Direction::new(x, T).unwrap(),
        T,
    )
    .unwrap();
    let circle = Circle::new(frame, radius, T).unwrap();
    make_edge(model, CircleCurve::new(circle).into(), (a, b), T)
        .unwrap()
        .shape
}

fn faces(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore(model, shape, Filter::OfType(ShapeType::Face)).unwrap()
}

fn edges(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore(model, shape, Filter::OfType(ShapeType::Edge)).unwrap()
}

fn surface(model: &Model, face: &Shape) -> SurfaceGeometry {
    let NodeData::Face(data) = model.node(face).unwrap().data() else {
        panic!("not a face");
    };
    model.geometry().surface(data.surface).unwrap().clone()
}

/// Points over the whole chart of a face's surface.
fn grid(surface: &SurfaceGeometry, steps: u32) -> Vec<Point> {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let mut out = Vec::new();
    for i in 0..=steps {
        for j in 0..=steps {
            let u = u0 + (u1 - u0) * f64::from(i) / f64::from(steps);
            let v = v0 + (v1 - v0) * f64::from(j) / f64::from(steps);
            out.push(surface.point_at(u, v, T).unwrap());
        }
    }
    out
}

/// A sheet of several faces checks clean but for its open border: the
/// faces share their inner edges, and exactly `open` edges bound it.
fn sound_sheet(model: &Model, shell: &Shape, open: usize) {
    let found = check(model, shell, T).unwrap();
    assert_eq!(found.problems.len(), 1, "{found:?}");
    let note = &found.problems[0];
    assert_eq!(note.severity, ogeom::algo::Severity::Suspect, "{found:?}");
    assert!(
        note.what
            .starts_with(&format!("{open} edge(s) used an odd number of times")),
        "{found:?}"
    );
}

/// The largest distance from points to a surface, by projection.
fn off(surface: &SurfaceGeometry, points: &[Point]) -> f64 {
    points
        .iter()
        .map(|p| project_on_surface(surface, *p, 24, T).unwrap().distance)
        .fold(0.0, f64::max)
}

#[test]
fn a_ruled_face_between_an_arc_and_a_line_is_bounded_by_them() {
    let mut model = Model::new();
    let a = arc(
        &mut model,
        Point::ORIGIN,
        Vector::new(0.0, 0.0, 1.0),
        Vector::new(1.0, 0.0, 0.0),
        5.0,
        0.0,
        PI,
    );
    let b = make_polygon(
        &mut model,
        &[Point::new(5.0, 0.0, 6.0), Point::new(-5.0, 0.0, 6.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let ruled = make_ruled(&mut model, &a, &b, T).unwrap().shape;
    assert_eq!(model.kind_of(&ruled).unwrap(), ShapeType::Face);
    assert!(check(&model, &ruled, T).unwrap().is_valid());

    // The long edges are the inputs themselves.
    let line = &edges(&model, &b)[0];
    let bounds = edges(&model, &ruled);
    assert_eq!(bounds.len(), 4);
    assert!(bounds.iter().any(|e| e.is_partner(&a)));
    assert!(bounds.iter().any(|e| e.is_partner(line)));

    // Every ruling is straight, from the arc at z = 0 to the line at
    // z = 6: the chart's v = 0 row on the circle, its v = 1 row on the
    // line, and the middle of each ruling on the line between its ends.
    let s = surface(&model, &faces(&model, &ruled)[0]);
    let ((u0, u1), (v0, v1)) = s.domain();
    for i in 0..=20 {
        let u = u0 + (u1 - u0) * f64::from(i) / 20.0;
        let p = s.point_at(u, v0, T).unwrap();
        let q = s.point_at(u, v1, T).unwrap();
        let m = s.point_at(u, 0.5 * (v0 + v1), T).unwrap();
        assert!(
            (p.x.hypot(p.y) - 5.0).abs() < 1e-9 && p.z.abs() < 1e-9,
            "{p:?}"
        );
        assert!(q.y.abs() < 1e-9 && (q.z - 6.0).abs() < 1e-9, "{q:?}");
        let along = (q - p) / (q - p).magnitude();
        let aside = (m - p) - along * (m - p).dot(along);
        assert!(
            aside.magnitude() < 1e-9,
            "ruling bends by {}",
            aside.magnitude()
        );
        assert!(m.z > 0.0 && m.z < 6.0, "{m:?}");
    }
}

#[test]
fn a_ruled_face_between_coplanar_segments_is_their_plane() {
    let mut model = Model::new();
    let a = make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(4.0, 0.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let b = make_polygon(
        &mut model,
        &[Point::new(1.0, 3.0, 0.0), Point::new(3.0, 3.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let ruled = make_ruled(&mut model, &a, &b, T).unwrap().shape;
    assert!(check(&model, &ruled, T).unwrap().is_valid());
    assert!(matches!(
        surface(&model, &faces(&model, &ruled)[0]),
        SurfaceGeometry::Plane(_)
    ));
    let props =
        ogeom::algo::surface_properties(&model, &ruled, ogeom::mesh::Deflection::default(), T)
            .unwrap();
    // The trapezoid: (4 + 2) / 2 * 3.
    assert!((props.mass - 9.0).abs() < 1e-9, "area {}", props.mass);
}

#[test]
fn a_ruled_face_between_skew_segments_is_the_bilinear_patch() {
    let mut model = Model::new();
    let a = make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(2.0, 0.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let b = make_polygon(
        &mut model,
        &[Point::new(0.0, 2.0, 1.0), Point::new(2.0, 2.0, -1.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let ruled = make_ruled(&mut model, &a, &b, T).unwrap().shape;
    assert!(check(&model, &ruled, T).unwrap().is_valid());
    // The hyperbolic paraboloid through the four corners: z = (1 - x) y / 2.
    let s = surface(&model, &faces(&model, &ruled)[0]);
    for p in grid(&s, 10) {
        assert!((p.z - (1.0 - p.x) * p.y / 2.0).abs() < 1e-12, "{p:?}");
    }
}

#[test]
fn a_ruled_face_between_circles_is_the_frustum() {
    let mut model = Model::new();
    let z = Vector::new(0.0, 0.0, 1.0);
    let x = Vector::new(1.0, 0.0, 0.0);
    let low = arc(&mut model, Point::ORIGIN, z, x, 2.0, 0.0, 2.0 * PI);
    let high = arc(
        &mut model,
        Point::new(0.0, 0.0, 3.0),
        z,
        x,
        1.0,
        0.0,
        2.0 * PI,
    );
    let low = make_wire(&mut model, &[low], T).unwrap().shape;
    let high = make_wire(&mut model, &[high], T).unwrap().shape;
    let ruled = make_ruled(&mut model, &low, &high, T).unwrap().shape;
    assert!(check(&model, &ruled, T).unwrap().is_valid());
    let s = surface(&model, &faces(&model, &ruled)[0]);
    for p in grid(&s, 16) {
        assert!((p.x.hypot(p.y) - (2.0 - p.z / 3.0)).abs() < 1e-9, "{p:?}");
    }
    // Lateral area of the frustum: pi (r1 + r2) slant.
    let props =
        ogeom::algo::surface_properties(&model, &ruled, ogeom::mesh::Deflection::default(), T)
            .unwrap();
    let want = PI * 3.0 * 10.0_f64.sqrt();
    assert!(
        (props.mass - want).abs() < want * 1e-6,
        "{} against {want}",
        props.mass
    );
}

#[test]
fn a_loft_through_three_open_arcs_passes_through_each() {
    let mut model = Model::new();
    let z = Vector::new(0.0, 0.0, 1.0);
    let x = Vector::new(1.0, 0.0, 0.0);
    let specs = [
        (0.0, 4.0, 0.0, PI / 2.0),
        (5.0, 6.0, 0.2, PI / 2.0 + 0.4),
        (10.0, 3.0, 0.0, 2.0),
    ];
    let sections: Vec<Shape> = specs
        .iter()
        .map(|&(h, r, a, b)| arc(&mut model, Point::new(0.0, 0.0, h), z, x, r, a, b))
        .collect();
    let loft = make_loft_surface(&mut model, &sections, false, &[], false, T)
        .unwrap()
        .shape;
    assert_eq!(model.kind_of(&loft).unwrap(), ShapeType::Face);
    assert!(check(&model, &loft, T).unwrap().is_valid());
    let s = surface(&model, &loft);
    for &(h, r, a, b) in &specs {
        let points: Vec<Point> = (0..=24)
            .map(|i| {
                let t = a + (b - a) * f64::from(i) / 24.0;
                Point::new(r * t.cos(), r * t.sin(), h)
            })
            .collect();
        let worst = off(&s, &points);
        assert!(
            worst < 1e-6,
            "section at z = {h} is {worst:.3e} off the loft"
        );
    }
    // The end sections bound it as the caller's own edges.
    let bounds = edges(&model, &loft);
    assert!(bounds.iter().any(|e| e.is_partner(&sections[0])));
    assert!(bounds.iter().any(|e| e.is_partner(&sections[2])));
}

#[test]
fn a_closed_loft_through_three_circles_comes_back_to_the_first() {
    let mut model = Model::new();
    let mut sections = Vec::new();
    let mut rings = Vec::new();
    for k in 0..3 {
        let th = 2.0 * PI * f64::from(k) / 3.0;
        let centre = Point::new(10.0 * th.cos(), 10.0 * th.sin(), 0.0);
        let normal = Vector::new(-th.sin(), th.cos(), 0.0);
        let radial = Vector::new(th.cos(), th.sin(), 0.0);
        let radius = [2.0, 3.0, 2.5][k as usize];
        let e = arc(&mut model, centre, normal, radial, radius, 0.0, 2.0 * PI);
        sections.push(make_wire(&mut model, &[e], T).unwrap().shape);
        rings.push((centre, normal, radial, radius));
    }
    let loft = make_loft_surface(&mut model, &sections, true, &[], false, T)
        .unwrap()
        .shape;
    assert_eq!(model.kind_of(&loft).unwrap(), ShapeType::Face);
    assert!(check(&model, &loft, T).unwrap().is_valid());
    let s = surface(&model, &loft);
    let ((u0, u1), (v0, v1)) = s.domain();
    for i in 0..=16 {
        let u = u0 + (u1 - u0) * f64::from(i) / 16.0;
        // The far end is the near end, slope and all.
        let (p, q) = (s.point_at(u, v0, T).unwrap(), s.point_at(u, v1, T).unwrap());
        assert!(p.distance(q) < 1e-12, "{p:?} against {q:?}");
        let (_, dp) = s.d1_at(u, v0, T).unwrap();
        let (_, dq) = s.d1_at(u, v1, T).unwrap();
        assert!(
            (dp - dq).magnitude() < 1e-9 * dp.magnitude(),
            "{dp:?} against {dq:?}"
        );
        // And it is the first circle.
        let (centre, normal, _, radius) = rings[0];
        assert!(((p - centre).magnitude() - radius).abs() < 1e-9);
        assert!((p - centre).dot(normal).abs() < 1e-9);
    }
    for &(centre, normal, radial, radius) in &rings {
        let side = normal.cross(radial);
        let points: Vec<Point> = (0..24)
            .map(|i| {
                let t = 2.0 * PI * f64::from(i) / 24.0;
                centre + (radial * t.cos() + side * t.sin()) * radius
            })
            .collect();
        let worst = off(&s, &points);
        assert!(worst < 1e-6, "a section is {worst:.3e} off the loop");
    }
}

#[test]
fn a_ruled_loft_has_a_face_per_span_sharing_the_middle_section() {
    let mut model = Model::new();
    let z = Vector::new(0.0, 0.0, 1.0);
    let x = Vector::new(1.0, 0.0, 0.0);
    let sections: Vec<Shape> = [(0.0, 3.0), (2.0, 4.0), (5.0, 2.0)]
        .iter()
        .map(|&(h, r)| arc(&mut model, Point::new(0.0, 0.0, h), z, x, r, 0.0, PI))
        .collect();
    let loft = make_loft_surface(&mut model, &sections, false, &[], true, T)
        .unwrap()
        .shape;
    assert_eq!(model.kind_of(&loft).unwrap(), ShapeType::Shell);
    sound_sheet(&model, &loft, 6);
    let fs = faces(&model, &loft);
    assert_eq!(fs.len(), 2);
    for f in &fs {
        assert!(edges(&model, f).iter().any(|e| e.is_partner(&sections[1])));
    }
}

#[test]
fn an_open_l_swept_along_a_quarter_circle_is_two_faces_and_no_caps() {
    let mut model = Model::new();
    // The quarter circle about (20, 0, 0) from the origin, heading -Y.
    let spine = arc(
        &mut model,
        Point::new(20.0, 0.0, 0.0),
        Vector::new(0.0, 0.0, 1.0),
        Vector::new(1.0, 0.0, 0.0),
        20.0,
        PI,
        1.5 * PI,
    );
    // An L in the plane square to the spine's start.
    let profile = make_polygon(
        &mut model,
        &[
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 2.0),
            Point::new(3.0, 0.0, 2.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let sheet = make_sweep_surface(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        T,
    )
    .unwrap()
    .shape;
    assert_eq!(model.kind_of(&sheet).unwrap(), ShapeType::Shell);
    sound_sheet(&model, &sheet, 6);
    let fs = faces(&model, &sheet);
    assert_eq!(fs.len(), 2, "one face per leg of the L, no caps");
    let legs = edges(&model, &profile);
    for leg in &legs {
        assert!(edges(&model, &sheet).iter().any(|e| e.is_partner(leg)));
    }
    // Along a planar circle the swept L turns about the circle's axis: in
    // the half-plane through the axis every point of the sheet lies on the
    // L, radius 19 for z in [0, 2] and z = 2 for radius in [17, 19].
    let mut widest: f64 = 0.0;
    for f in &fs {
        assert!(!matches!(surface(&model, f), SurfaceGeometry::Plane(_)));
        for p in grid(&surface(&model, f), 24) {
            let rho = (p.x - 20.0).hypot(p.y);
            let upright = (rho - 19.0).abs().max((p.z.clamp(0.0, 2.0) - p.z).abs());
            let flat = (p.z - 2.0).abs().max((rho.clamp(17.0, 19.0) - rho).abs());
            let d = upright.min(flat);
            assert!(d < 1e-6, "{p:?} is {d:.3e} off the swept L");
            widest = widest.max(p.y.atan2(20.0 - p.x).abs());
        }
    }
    // A quarter turn.
    assert!((widest - PI / 2.0).abs() < 1e-9, "turns {widest}");
    // Pappus: the upright leg sweeps 2 x 19 x pi / 2, the flat one a
    // quarter of the annulus between radii 17 and 19.
    let area =
        ogeom::algo::surface_properties(&model, &sheet, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    let want = 19.0 * PI + 18.0 * PI;
    assert!((area - want).abs() < want * 1e-6, "{area} against {want}");
}

#[test]
fn a_profile_between_two_rails_widens_with_them() {
    let mut model = Model::new();
    let rail_a = make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(0.0, 10.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let rail_b = make_polygon(
        &mut model,
        &[Point::new(2.0, 0.0, 0.0), Point::new(4.0, 10.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    // A half circle hanging below the rails' starts.
    let profile = arc(
        &mut model,
        Point::new(1.0, 0.0, 0.0),
        Vector::new(0.0, -1.0, 0.0),
        Vector::new(-1.0, 0.0, 0.0),
        1.0,
        0.0,
        PI,
    );
    let sheet = make_sweep_two_rails(&mut model, &profile, &rail_a, &rail_b, T)
        .unwrap()
        .shape;
    assert_eq!(model.kind_of(&sheet).unwrap(), ShapeType::Face);
    assert!(check(&model, &sheet, T).unwrap().is_valid());
    assert!(edges(&model, &sheet).iter().any(|e| e.is_partner(&profile)));
    // At height y the rails stand at x = 0 and x = 2 + y / 5, both a
    // fraction y / 10 along: the section there is the half circle on that
    // diameter, upright.
    let s = surface(&model, &sheet);
    for p in grid(&s, 20) {
        let (left, right) = (0.0, 2.0 + p.y / 5.0);
        let centre = Point::new(f64::midpoint(left, right), p.y, 0.0);
        let d = (p.distance(centre) - (right - left) / 2.0).abs();
        assert!(d < 1e-6, "{p:?} is {d:.3e} off the widening half circle");
        assert!(p.z < 1e-9, "{p:?} above the rails' plane");
    }
}

/// The arc of the circle through `a`, `b` and `c`, from `a` through `b` to
/// `c`, with the circle itself and the arc's end angle (it starts at 0).
fn arc_through(model: &mut Model, a: Point, b: Point, c: Point) -> (Shape, Circle, f64) {
    let circle = Circle::through(a, b, c, T).unwrap();
    let end = angle_on(&circle, c);
    assert!(angle_on(&circle, b) < end);
    let edge = make_edge(model, CircleCurve::new(circle).into(), (0.0, end), T)
        .unwrap()
        .shape;
    (edge, circle, end)
}

/// Points along an arc of `circle` from angle `a` to `b`.
fn on_circle(circle: &Circle, a: f64, b: f64, count: u32) -> Vec<Point> {
    let f = circle.frame();
    (0..=count)
        .map(|i| {
            let t = a + (b - a) * f64::from(i) / f64::from(count);
            f.origin() + (f.x().vector() * t.cos() + f.y().vector() * t.sin()) * circle.radius()
        })
        .collect()
}

/// The angle of `p` about `circle`'s centre, from its `x`, in `[0, 2 pi)`.
fn angle_on(circle: &Circle, p: Point) -> f64 {
    let f = circle.frame();
    let d = p - f.origin();
    let t = d.dot(f.y().vector()).atan2(d.dot(f.x().vector()));
    if t < 0.0 { t + 2.0 * PI } else { t }
}

/// The distance from `p` to a circle.
fn off_circle(circle: &Circle, p: Point) -> f64 {
    let f = circle.frame();
    let d = p - f.origin();
    let up = d.dot(f.z().vector());
    let flat = (d - f.z().vector() * up).magnitude();
    (flat - circle.radius()).hypot(up)
}

#[test]
fn a_loft_through_three_open_arcs_follows_two_guides_along_their_ends() {
    let mut model = Model::new();
    let z = Vector::new(0.0, 0.0, 1.0);
    let x = Vector::new(1.0, 0.0, 0.0);
    let specs = [
        (0.0, 4.0, 0.0, PI / 2.0),
        (5.0, 6.0, 0.2, PI / 2.0 + 0.4),
        (10.0, 3.0, 0.0, 2.0),
    ];
    let sections: Vec<Shape> = specs
        .iter()
        .map(|&(h, r, a, b)| arc(&mut model, Point::new(0.0, 0.0, h), z, x, r, a, b))
        .collect();
    let at = |r: f64, t: f64, h: f64| Point::new(r * t.cos(), r * t.sin(), h);
    // A cubic spline through the arcs' starts, bowing out between them,
    // and a circular arc through their ends.
    let starts = [
        at(4.0, 0.0, 0.0),
        Point::new(6.5, -0.5, 2.5),
        at(6.0, 0.2, 5.0),
        Point::new(6.0, -0.5, 7.5),
        at(3.0, 0.0, 10.0),
    ];
    let spline = interpolate(&starts, 3, Spacing::Centripetal, T).unwrap();
    let domain = spline.knots().domain();
    let first = make_edge(&mut model, Curve::BSpline(spline.clone()), domain, T)
        .unwrap()
        .shape;
    let ends: Vec<Point> = specs.iter().map(|&(h, r, _, b)| at(r, b, h)).collect();
    let (second, circle, end) = arc_through(&mut model, ends[0], ends[1], ends[2]);
    let unguided = make_loft_surface(&mut model, &sections, false, &[], false, T)
        .unwrap()
        .shape;
    let loft = make_loft_surface(&mut model, &sections, false, &[first, second], false, T)
        .unwrap()
        .shape;
    assert_eq!(model.kind_of(&loft).unwrap(), ShapeType::Face);
    assert!(check(&model, &loft, T).unwrap().is_valid());
    let s = surface(&model, &loft);
    for &(h, r, a, b) in &specs {
        let points: Vec<Point> = (0..=24)
            .map(|i| at(r, a + (b - a) * f64::from(i) / 24.0, h))
            .collect();
        let worst = off(&s, &points);
        assert!(
            worst < 1e-6,
            "section at z = {h} is {worst:.3e} off the loft"
        );
    }
    // Each guide lies on the loft within the stated ten confusions, and
    // the unguided loft is off them: well off the bowed spline, and off
    // the arc by far more than that tolerance.
    let guide_points = [
        (0..=48)
            .map(|i| {
                let t = domain.0 + (domain.1 - domain.0) * f64::from(i) / 48.0;
                spline.point_at(t, T).unwrap()
            })
            .collect::<Vec<Point>>(),
        on_circle(&circle, 0.0, end, 48),
    ];
    let plain = surface(&model, &unguided);
    for (points, away) in guide_points.iter().zip([0.3, 1e-3]) {
        let worst = off(&s, points);
        assert!(worst < 1e-6, "a guide is {worst:.3e} off the loft");
        let apart = off(&plain, points);
        assert!(
            apart > away,
            "the unguided loft is already {apart:.3e} from a guide"
        );
    }
    // The end sections still bound it as the caller's own edges.
    let bounds = edges(&model, &loft);
    assert!(bounds.iter().any(|e| e.is_partner(&sections[0])));
    assert!(bounds.iter().any(|e| e.is_partner(&sections[2])));
}

#[test]
fn a_loft_between_two_lines_bulges_to_follow_a_curved_guide() {
    let mut model = Model::new();
    let line = |model: &mut Model, h: f64| {
        make_polygon(
            model,
            &[Point::new(-5.0, 0.0, h), Point::new(5.0, 0.0, h)],
            false,
            T,
        )
        .unwrap()
        .shape
    };
    let lower = line(&mut model, 0.0);
    let upper = line(&mut model, 10.0);
    let (guide, circle, end) = arc_through(
        &mut model,
        Point::new(0.0, 0.0, 0.0),
        Point::new(0.0, 3.0, 5.0),
        Point::new(0.0, 0.0, 10.0),
    );
    // The circle through those three points: centre (0, -8/3, 5), radius
    // 17/3.
    assert!((circle.radius() - 17.0 / 3.0).abs() < 1e-12);
    let loft = make_loft_surface(&mut model, &[lower, upper], false, &[guide], false, T)
        .unwrap()
        .shape;
    assert_eq!(model.kind_of(&loft).unwrap(), ShapeType::Face);
    assert!(check(&model, &loft, T).unwrap().is_valid());
    let s = surface(&model, &loft);
    let along = on_circle(&circle, 0.0, end, 64);
    let worst = off(&s, &along);
    assert!(worst < 1e-6, "the guide is {worst:.3e} off the loft");
    for h in [0.0, 10.0] {
        let points: Vec<Point> = (0..=20)
            .map(|i| Point::new(-5.0 + f64::from(i) * 0.5, 0.0, h))
            .collect();
        let worst = off(&s, &points);
        assert!(
            worst < 1e-6,
            "the line at z = {h} is {worst:.3e} off the loft"
        );
    }
    // The bulge: out of the lines' plane by up to the guide's three, and
    // no farther.
    let reach = grid(&s, 32).iter().map(|p| p.y).fold(0.0, f64::max);
    assert!(
        reach > 2.9 && reach < 3.0 + 1e-9,
        "the loft reaches y = {reach}"
    );
}

#[test]
fn a_guide_crossing_lines_off_centre_keeps_them_as_the_bounding_edges() {
    let mut model = Model::new();
    let line = |model: &mut Model, h: f64| {
        make_polygon(
            model,
            &[Point::new(-5.0, 0.0, h), Point::new(5.0, 0.0, h)],
            false,
            T,
        )
        .unwrap()
        .shape
    };
    let lower = line(&mut model, 0.0);
    let upper = line(&mut model, 10.0);
    // Crossing the lower line at 0.4 of its length and the upper at 0.7:
    // each line is paced anew so the guide crosses both at one parameter.
    let (guide, circle, end) = arc_through(
        &mut model,
        Point::new(-1.0, 0.0, 0.0),
        Point::new(0.5, 2.0, 5.0),
        Point::new(2.0, 0.0, 10.0),
    );
    let loft = make_loft_surface(
        &mut model,
        &[lower.clone(), upper.clone()],
        false,
        &[guide],
        false,
        T,
    )
    .unwrap()
    .shape;
    assert!(check(&model, &loft, T).unwrap().is_valid());
    let s = surface(&model, &loft);
    let worst = off(&s, &on_circle(&circle, 0.0, end, 64));
    assert!(worst < 1e-6, "the guide is {worst:.3e} off the loft");
    let mut bounding = 0;
    let bounds = edges(&model, &loft);
    for (shape, h) in [(&lower, 0.0), (&upper, 10.0)] {
        let points: Vec<Point> = (0..=20)
            .map(|i| Point::new(-5.0 + f64::from(i) * 0.5, 0.0, h))
            .collect();
        let worst = off(&s, &points);
        assert!(
            worst < 1e-6,
            "the line at z = {h} is {worst:.3e} off the loft"
        );
        let own = &edges(&model, shape)[0];
        bounding += usize::from(bounds.iter().any(|e| e.is_partner(own)));
        // Its image in the loft's chart lands where the line is at the
        // same parameter.
        let gap = pcurve_gap(&model, &loft, own);
        assert!(gap < 1e-6, "the line's image strays {gap:.3e} from it");
    }
    assert_eq!(bounding, 2, "the lines bound the loft as their own edges");
}

/// The largest distance between an edge's curve and the face's surface
/// at the edge's image in the face's chart, at the same parameter.
fn pcurve_gap(model: &Model, face: &Shape, edge: &Shape) -> f64 {
    let id = model.node(face).unwrap().data().as_face().unwrap().surface;
    let s = model.geometry().surface(id).unwrap();
    let data = model.node(edge).unwrap().data().as_edge().unwrap();
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        panic!("an edge with no curve");
    };
    let Some(EdgeRepr::PCurve { curve: image, .. }) = data.pcurve_on(id) else {
        panic!("an edge with no image on the face");
    };
    let curve = model.geometry().curve(*curve).unwrap();
    let image = model.geometry().pcurve(*image).unwrap();
    (0..=200)
        .map(|i| {
            let t = range.0 + (range.1 - range.0) * f64::from(i) / 200.0;
            let q = image.point_at(t, T).unwrap();
            s.point_at(q.x, q.y, T)
                .unwrap()
                .distance(curve.point_at(t, T).unwrap())
        })
        .fold(0.0, f64::max)
}

#[test]
fn a_guided_loft_through_circles_closes_on_its_seam_guide() {
    let mut model = Model::new();
    let z = Vector::new(0.0, 0.0, 1.0);
    let x = Vector::new(1.0, 0.0, 0.0);
    let rings = [(0.0, 3.0), (5.0, 4.0), (10.0, 3.0)];
    let sections: Vec<Shape> = rings
        .iter()
        .map(|&(h, r)| arc(&mut model, Point::new(0.0, 0.0, h), z, x, r, 0.0, 2.0 * PI))
        .collect();
    let at = |r: f64, t: f64, h: f64| Point::new(r * t.cos(), r * t.sin(), h);
    // The seam guide crosses the circles away from their starts, at
    // different angles, and the second guide at angles that are not one
    // fraction of the turn past the first: sections and guides both take
    // a new pace.
    let pick = |angles: [f64; 3]| -> Vec<Point> {
        rings
            .iter()
            .zip(angles)
            .map(|(&(h, r), t)| at(r, t, h))
            .collect()
    };
    let p = pick([0.3, 0.5, 0.4]);
    let (seam, seam_circle, seam_end) = arc_through(&mut model, p[0], p[1], p[2]);
    let q = pick([2.0, 2.6, 2.2]);
    let (other, other_circle, other_end) = arc_through(&mut model, q[0], q[1], q[2]);
    let loft = make_loft_surface(&mut model, &sections, false, &[seam, other], false, T)
        .unwrap()
        .shape;
    assert_eq!(model.kind_of(&loft).unwrap(), ShapeType::Face);
    assert!(check(&model, &loft, T).unwrap().is_valid());
    let s = surface(&model, &loft);
    for &(h, r) in &rings {
        let points: Vec<Point> = (0..48)
            .map(|i| at(r, 2.0 * PI * f64::from(i) / 48.0, h))
            .collect();
        let worst = off(&s, &points);
        assert!(
            worst < 1e-6,
            "the circle at z = {h} is {worst:.3e} off the loft"
        );
    }
    for (circle, end) in [(seam_circle, seam_end), (other_circle, other_end)] {
        let worst = off(&s, &on_circle(&circle, 0.0, end, 48));
        assert!(worst < 1e-6, "a guide is {worst:.3e} off the loft");
    }
    // The chart's two ends along the sections are one curve, on the seam
    // guide.
    let ((u0, u1), (v0, v1)) = s.domain();
    for j in 0..=16 {
        let v = v0 + (v1 - v0) * f64::from(j) / 16.0;
        let (a, b) = (s.point_at(u0, v, T).unwrap(), s.point_at(u1, v, T).unwrap());
        assert!(a.distance(b) < 1e-9, "{a:?} against {b:?}");
        let gap = off_circle(&seam_circle, a);
        assert!(gap < 1e-6, "the seam is {gap:.3e} off its guide");
    }
}

/// A segment along `x` from 0 to 1 at height `h`.
fn rung(model: &mut Model, h: f64) -> Shape {
    make_polygon(
        model,
        &[Point::new(0.0, 0.0, h), Point::new(1.0, 0.0, h)],
        false,
        T,
    )
    .unwrap()
    .shape
}

/// A segment from `x0` at height 0 to `x1` at height 2.
fn climb(model: &mut Model, x0: f64, x1: f64) -> Shape {
    make_polygon(
        model,
        &[Point::new(x0, 0.0, 0.0), Point::new(x1, 0.0, 2.0)],
        false,
        T,
    )
    .unwrap()
    .shape
}

#[test]
fn a_loft_surface_refuses_a_guide_that_misses_a_section_by_name() {
    let mut model = Model::new();
    let sections = [
        rung(&mut model, 0.0),
        rung(&mut model, 1.0),
        rung(&mut model, 2.0),
    ];
    // Up the sections' starts, stopping half a unit short of the last.
    let guide = make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(0.0, 0.0, 1.5)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let err = make_loft_surface(&mut model, &sections, false, &[guide], false, T).unwrap_err();
    assert!(
        err.to_string()
            .contains("guide 0 misses section 2 by 5.000e-1"),
        "{err}"
    );
}

#[test]
fn a_guided_loft_surface_refuses_what_it_does_not_build_by_name() {
    let mut model = Model::new();
    let three = [
        rung(&mut model, 0.0),
        rung(&mut model, 1.0),
        rung(&mut model, 2.0),
    ];
    let straight = [climb(&mut model, 0.5, 0.5)];
    let refused = |model: &mut Model, sections: &[Shape], closed, guides: &[Shape], ruled| {
        make_loft_surface(model, sections, closed, guides, ruled, T)
            .unwrap_err()
            .to_string()
    };
    let err = refused(&mut model, &three, false, &straight, true);
    assert!(err.contains("follows no guide curves"), "{err}");
    let err = refused(&mut model, &three, true, &straight, false);
    assert!(
        err.contains("a closed loft surface does not follow guide curves"),
        "{err}"
    );
    let bent = make_polygon(
        &mut model,
        &[
            Point::new(0.0, 0.0, 2.0),
            Point::new(0.5, 0.0, 2.0),
            Point::new(1.0, 0.0, 2.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let err = refused(
        &mut model,
        &[three[0].clone(), bent],
        false,
        &straight,
        false,
    );
    assert!(
        err.contains("takes sections of one edge each; section 1 has 2"),
        "{err}"
    );
    // Two guides crossing each other between the sections.
    let crossed = [climb(&mut model, 0.2, 0.8), climb(&mut model, 0.8, 0.2)];
    let err = refused(&mut model, &three, false, &crossed, false);
    assert!(
        err.contains(
            "guides 0 and 1 cross section 1 in another order than section 0, or at one point"
        ),
        "{err}"
    );
    // A guide that crosses the first section at its start and the others
    // inside them.
    let slanted = [climb(&mut model, 0.0, 1.0)];
    let err = refused(&mut model, &three, false, &slanted, false);
    assert!(
        err.contains("guide 0 crosses some sections at an end and others inside them"),
        "{err}"
    );
}

#[test]
fn a_loft_surface_refuses_sections_of_different_edge_counts() {
    let mut model = Model::new();
    let a = make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(1.0, 0.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let b = make_polygon(
        &mut model,
        &[
            Point::new(0.0, 0.0, 1.0),
            Point::new(0.5, 0.5, 1.0),
            Point::new(1.0, 0.0, 1.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let err = make_ruled(&mut model, &a, &b, T).unwrap_err();
    assert!(err.to_string().contains("pair edge for edge"), "{err}");
}

#[test]
fn a_sweep_surface_refuses_a_closed_spine_by_name() {
    let mut model = Model::new();
    let spine = arc(
        &mut model,
        Point::ORIGIN,
        Vector::new(0.0, 0.0, 1.0),
        Vector::new(1.0, 0.0, 0.0),
        10.0,
        0.0,
        2.0 * PI,
    );
    let profile = make_polygon(
        &mut model,
        &[Point::new(11.0, 0.0, 0.0), Point::new(11.0, 0.0, 1.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let err = make_sweep_surface(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        T,
    )
    .unwrap_err();
    assert!(
        err.to_string().contains("closed spine is not built"),
        "{err}"
    );
}

#[test]
fn a_two_rail_sweep_refuses_a_profile_off_the_rails() {
    let mut model = Model::new();
    let rail_a = make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(0.0, 10.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let rail_b = make_polygon(
        &mut model,
        &[Point::new(2.0, 0.0, 0.0), Point::new(2.0, 10.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let profile = make_polygon(
        &mut model,
        &[Point::new(0.0, 0.0, 1.0), Point::new(2.0, 0.0, 1.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let err = make_sweep_two_rails(&mut model, &profile, &rail_a, &rail_b, T).unwrap_err();
    assert!(
        err.to_string().contains("must sit on the rails' starts"),
        "{err}"
    );
}

/// A closed planar polygon at height `z`.
fn polygon(model: &mut Model, corners: &[(f64, f64)], z: f64) -> Shape {
    let points: Vec<Point> = corners.iter().map(|&(x, y)| Point::new(x, y, z)).collect();
    make_polygon(model, &points, true, T).unwrap().shape
}

#[test]
fn a_ruled_sheet_between_two_squares_sews_with_their_caps_into_a_closed_shell() {
    let mut model = Model::new();
    let low = polygon(
        &mut model,
        &[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)],
        0.0,
    );
    let high = polygon(
        &mut model,
        &[(0.5, 0.0), (2.5, 0.5), (2.0, 2.5), (0.0, 2.0)],
        3.0,
    );
    let sheet = make_ruled(&mut model, &low, &high, T).unwrap().shape;
    // Four walls, their eight section edges open.
    sound_sheet(&model, &sheet, 8);
    let mut all = faces(&model, &sheet);
    for (wire, z, up) in [(&low, 0.0, -1.0), (&high, 3.0, 1.0)] {
        let plane = ogeom::math::Plane::new(
            Frame::new(
                Point::new(0.0, 0.0, z),
                Direction::new(Vector::new(0.0, 0.0, up), T).unwrap(),
                Direction::X,
                T,
            )
            .unwrap(),
        );
        let cap = ogeom::algo::make_face(
            &mut model,
            ogeom::geom::PlaneSurface::new(plane).into(),
            std::slice::from_ref(wire),
            T,
        )
        .unwrap()
        .shape;
        all.push(cap);
    }
    let sewn = ogeom::algo::sew(&mut model, &all, T).unwrap();
    assert_eq!(sewn.shells.len(), 1);
    assert!(ogeom::algo::is_shell_closed(&model, &sewn.shells[0]).unwrap());
}

#[test]
fn a_smooth_loft_through_closed_squares_is_a_face_per_side() {
    let mut model = Model::new();
    let square = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
    let sections: Vec<Shape> = [(0.0, 1.0), (2.0, 2.0), (4.0, 1.5)]
        .iter()
        .map(|&(z, k)| {
            let corners: Vec<(f64, f64)> = square.iter().map(|&(x, y)| (x * k, y * k)).collect();
            polygon(&mut model, &corners, z)
        })
        .collect();
    let loft = make_loft_surface(&mut model, &sections, false, &[], false, T)
        .unwrap()
        .shape;
    sound_sheet(&model, &loft, 8);
    // Every section lies on the sheet: each point of a square's side on
    // one of the four faces, to the projection's own precision.
    let surfaces: Vec<SurfaceGeometry> = faces(&model, &loft)
        .iter()
        .map(|f| surface(&model, f))
        .collect();
    for (z, k) in [(0.0, 1.0), (2.0, 2.0), (4.0, 1.5)] {
        for i in 0..4 {
            let (a, b) = (square[i], square[(i + 1) % 4]);
            for j in 0..=8 {
                let f = f64::from(j) / 8.0;
                let p = Point::new(k * (a.0 + (b.0 - a.0) * f), k * (a.1 + (b.1 - a.1) * f), z);
                let d = surfaces
                    .iter()
                    .map(|s| off(s, &[p]))
                    .fold(f64::INFINITY, f64::min);
                assert!(d < 1e-8, "{p:?} is {d:.3e} off the sheet");
            }
        }
    }
}

#[test]
fn a_closed_loft_through_open_arcs_is_one_band_round_the_loop() {
    let mut model = Model::new();
    let mut arcs = Vec::new();
    for k in 0..4 {
        let th = 2.0 * PI * f64::from(k) / 4.0;
        let centre = Point::new(10.0 * th.cos(), 10.0 * th.sin(), 0.0);
        let normal = Vector::new(-th.sin(), th.cos(), 0.0);
        let radial = Vector::new(th.cos(), th.sin(), 0.0);
        arcs.push(arc(&mut model, centre, normal, radial, 2.0, 0.0, PI / 2.0));
    }
    let band = make_loft_surface(&mut model, &arcs, true, &[], false, T)
        .unwrap()
        .shape;
    assert_eq!(model.kind_of(&band).unwrap(), ShapeType::Face);
    let found = check(&model, &band, T).unwrap();
    assert!(found.is_valid(), "{found:?}");
    let s = surface(&model, &band);
    for e in &arcs {
        let (curve, range) = {
            let NodeData::Edge(data) = model.node(e).unwrap().data() else {
                panic!("not an edge");
            };
            let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                panic!("no curve");
            };
            (model.geometry().curve(*curve).unwrap().clone(), *range)
        };
        let points: Vec<Point> = (0..=16)
            .map(|i| {
                use ogeom::geom::Curve3d as _;
                curve
                    .point_at(range.0 + (range.1 - range.0) * f64::from(i) / 16.0, T)
                    .unwrap()
            })
            .collect();
        let worst = off(&s, &points);
        assert!(worst < 1e-6, "an arc is {worst:.3e} off the band");
    }
}

#[test]
fn a_fixed_sweep_along_a_line_is_the_extrusion() {
    let mut model = Model::new();
    let spine = make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(0.0, 3.0, 4.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let profile = arc(
        &mut model,
        Point::ORIGIN,
        Vector::new(0.0, 0.0, 1.0),
        Vector::new(1.0, 0.0, 0.0),
        2.0,
        0.0,
        PI,
    );
    let sheet = make_sweep_surface(&mut model, &profile, &spine, &PipeLaw::Fixed, T)
        .unwrap()
        .shape;
    assert!(check(&model, &sheet, T).unwrap().is_valid());
    // Every point is the arc moved along (0, 3, 4): back along it to
    // z = 0, it lies on the circle of radius 2.
    let s = surface(&model, &sheet);
    for p in grid(&s, 12) {
        let f = p.z / 4.0;
        let q = Point::new(p.x, p.y - 3.0 * f, 0.0);
        assert!((q.x.hypot(q.y) - 2.0).abs() < 1e-9, "{p:?}");
        assert!(q.y > -1e-9 && (-1e-9..=4.0 + 1e-9).contains(&p.z));
    }
}
