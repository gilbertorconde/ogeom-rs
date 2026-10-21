//! Pipe shells whose section scales and twists along the spine, measured
//! against closed forms: a linear scale down a straight spine is a
//! frustum or a pyramid, exactly; a twist keeps every slice the profile
//! turned, so the volume is the slices' areas summed along the spine.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_edge, make_face, make_polygon, make_wire, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve, Curve3d as _, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::offset::{PipeEvolution, PipeLaw, make_pipe_shell_evolving};
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// The skin's tolerance for the skinned sweeps.
const SKIN: f64 = 1e-3;

fn line(model: &mut Model, a: Point, b: Point) -> Shape {
    let curve = LineCurve::segment(a, b, T).unwrap();
    let range = curve.domain();
    make_edge(model, Curve::Line(curve), range, T)
        .unwrap()
        .shape
}

/// A circle of `radius` about `centre` in the `xy` plane.
fn circle_wire(model: &mut Model, centre: Point, radius: f64) -> Shape {
    let frame = Frame::new(centre, Direction::Z, Direction::X, T).unwrap();
    let circle = Circle::new(frame, radius, T).unwrap();
    let edge = make_edge(model, CircleCurve::new(circle).into(), (0.0, 2.0 * PI), T)
        .unwrap()
        .shape;
    make_wire(model, &[edge], T).unwrap().shape
}

fn square_wire(model: &mut Model, half: f64) -> Shape {
    let pts = [(-half, -half), (half, -half), (half, half), (-half, half)]
        .map(|(x, y)| Point::new(x, y, 0.0));
    make_polygon(model, &pts, true, T).unwrap().shape
}

fn xy_face(model: &mut Model, wires: &[Shape]) -> Shape {
    let plane = PlaneSurface::new(Plane::new(Frame::WORLD));
    make_face(model, plane.into(), wires, T).unwrap().shape
}

/// Volume, and whether it was integrated exactly.
fn measure(model: &Model, shape: &Shape) -> (f64, bool) {
    let p = volume_properties(model, shape, Deflection::default(), T).unwrap();
    (p.mass, p.deflection == 0.0)
}

fn surfaces(model: &Model, shape: &Shape) -> Vec<SurfaceGeometry> {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .map(|f| {
            let Some(NodeData::Face(d)) = model.node(f).map(|n| n.data()) else {
                panic!("a face holds no face data");
            };
            model.geometry().surface(d.surface).unwrap().clone()
        })
        .collect()
}

/// A section scaled linearly from one to `s` over a run `h` has area
/// `A (1 + (s - 1) λ)²`, so the solid holds `A h (1 + s + s²) / 3`
/// whatever the section turns through, wherever its centroid stays on a
/// spine it stands square to.
fn scaled_volume(area: f64, h: f64, s: f64) -> f64 {
    area * h * (1.0 + s + s * s) / 3.0
}

fn evolve(scale: f64, twist: f64) -> PipeEvolution {
    PipeEvolution { scale, twist }
}

#[test]
fn a_square_scaled_down_a_line_is_the_pyramid_s_frustum_exactly() {
    let (half, h, s) = (1.5, 6.0, 0.4);
    let mut model = Model::new();
    let wire = square_wire(&mut model, half);
    let profile = xy_face(&mut model, &[wire]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
    let solid = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        evolve(s, 0.0),
        SKIN,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let all = surfaces(&model, &solid);
    assert_eq!(all.len(), 6);
    assert!(all.iter().all(|s| matches!(s, SurfaceGeometry::Plane(_))));
    let (v, exact) = measure(&model, &solid);
    assert!(exact, "a solid of planes integrates exactly");
    let expected = scaled_volume(4.0 * half * half, h, s);
    assert!(
        (v - expected).abs() < expected * 1e-12,
        "{v} against {expected}"
    );
}

#[test]
fn a_disc_scaled_down_a_line_is_the_cone_s_frustum_and_to_naught_the_cone() {
    let (r, h) = (2.0, 5.0);
    for s in [2.5, 0.0] {
        let mut model = Model::new();
        let wire = circle_wire(&mut model, Point::ORIGIN, r);
        let profile = xy_face(&mut model, &[wire]);
        let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
        let solid = make_pipe_shell_evolving(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            evolve(s, 0.0),
            SKIN,
            T,
        )
        .unwrap()
        .shape;
        let found = check(&model, &solid, T).unwrap();
        assert!(found.is_valid(), "scale {s}: {found}");
        assert!(
            surfaces(&model, &solid)
                .iter()
                .any(|f| matches!(f, SurfaceGeometry::Cone(_))),
            "scale {s}: the wall is the cone"
        );
        let (v, _) = measure(&model, &solid);
        let expected = scaled_volume(PI * r * r, h, s);
        assert!(
            (v - expected).abs() < expected * 1e-9,
            "scale {s}: {v} against {expected}"
        );
    }
}

#[test]
fn a_square_scaled_to_naught_is_the_pyramid() {
    let (half, h) = (1.0, 3.0);
    let mut model = Model::new();
    let wire = square_wire(&mut model, half);
    let profile = xy_face(&mut model, &[wire]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
    let solid = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        evolve(0.0, 0.0),
        SKIN,
        T,
    )
    .unwrap()
    .shape;
    assert!(check(&model, &solid, T).unwrap().is_valid());
    assert_eq!(surfaces(&model, &solid).len(), 5, "four sides and a base");
    let (v, _) = measure(&model, &solid);
    let expected = 4.0 * half * half * h / 3.0;
    assert!(
        (v - expected).abs() < expected * 1e-12,
        "{v} against {expected}"
    );
}

/// The chord deflection a skinned sweep's volume is measured to.
const CHORD: f64 = 1e-4;

/// A skinned sweep's volume, integrated over its mesh at [`CHORD`].
fn measure_skinned(model: &Model, shape: &Shape) -> f64 {
    let fine = Deflection {
        chord: CHORD,
        ..Deflection::default()
    };
    volume_properties(model, shape, fine, T).unwrap().mass
}

/// The measured volume of a skinned sweep strays from the swept solid's by
/// at most the skin's tolerance and the mesh's chord over the walls' area:
/// their sum times a bound on the area, the perimeter at its largest times
/// the longest path a profile point runs.
fn skin_bound(perimeter: f64, longest_path: f64) -> f64 {
    (SKIN + CHORD) * perimeter * longest_path
}

#[test]
fn a_twisted_square_down_a_line_holds_the_prism_s_volume() {
    let (half, h, twist) = (1.0, 4.0, PI / 2.0);
    for s in [1.0, 0.5] {
        let mut model = Model::new();
        let wire = square_wire(&mut model, half);
        let profile = xy_face(&mut model, &[wire]);
        let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
        let solid = make_pipe_shell_evolving(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            evolve(s, twist),
            SKIN,
            T,
        )
        .unwrap()
        .shape;
        let found = check(&model, &solid, T).unwrap();
        assert!(found.is_usable(), "scale {s}: {found}");
        // The far end is the square turned a quarter turn and scaled: its
        // corners stand on the diagonals again.
        let corner = Point::new(half * s, half * s, h);
        let vertices = explore_unique(&model, &solid, ShapeType::Vertex).unwrap();
        let nearest = vertices
            .iter()
            .map(|v| {
                let Some(NodeData::Vertex(d)) = model.node(v).map(|n| n.data()) else {
                    panic!("a vertex holds no point");
                };
                d.point.distance(corner)
            })
            .fold(f64::INFINITY, f64::min);
        assert!(
            nearest < 10.0 * T.confusion(),
            "scale {s}: corner off by {nearest}"
        );
        let v = measure_skinned(&model, &solid);
        let expected = scaled_volume(4.0 * half * half, h, s);
        let radius = half * 2.0_f64.sqrt();
        let bound = skin_bound(8.0 * half, h.hypot(twist * radius));
        assert!(
            (v - expected).abs() < bound,
            "scale {s}: {v} against {expected}, bound {bound}"
        );
    }
}

/// A square centred on a quarter circle, scaled as it goes: its centroid
/// stays on the spine, so by Pappus the solid holds the sections' areas
/// summed along the spine's length.
#[test]
fn a_square_scaled_round_an_arc_holds_pappus_s_volume() {
    let (big, half, s) = (10.0, 1.0, 0.5);
    let mut model = Model::new();
    // The arc in the xz plane from the origin, leaving it along z.
    let frame = Frame::new(
        Point::new(big, 0.0, 0.0),
        Direction::new(Vector::new(0.0, 1.0, 0.0), T).unwrap(),
        Direction::new(Vector::new(-1.0, 0.0, 0.0), T).unwrap(),
        T,
    )
    .unwrap();
    let arc = Circle::new(frame, big, T).unwrap();
    let spine = make_edge(&mut model, CircleCurve::new(arc).into(), (0.0, PI / 2.0), T)
        .unwrap()
        .shape;
    let wire = square_wire(&mut model, half);
    let profile = xy_face(&mut model, &[wire]);
    let solid = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        evolve(s, 0.0),
        SKIN,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_usable(), "{found}");
    let length = big * PI / 2.0;
    let v = measure_skinned(&model, &solid);
    let expected = scaled_volume(4.0 * half * half, length, s);
    let reach = half * 2.0_f64.sqrt();
    let bound = skin_bound(8.0 * half, length * (1.0 + reach / big));
    assert!(
        (v - expected).abs() < bound,
        "{v} against {expected}, bound {bound}"
    );
}

#[test]
fn a_holed_square_scaled_down_a_line_keeps_its_hole() {
    let (half, r, h, s) = (2.0, 0.75, 5.0, 0.5);
    let mut model = Model::new();
    let outer = square_wire(&mut model, half);
    let hole = circle_wire(&mut model, Point::ORIGIN, r);
    let profile = xy_face(&mut model, &[outer, hole]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
    let solid = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        evolve(s, 0.0),
        SKIN,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let (v, _) = measure(&model, &solid);
    let expected = scaled_volume(4.0 * half * half - PI * r * r, h, s);
    assert!(
        (v - expected).abs() < expected * 1e-9,
        "{v} against {expected}"
    );
}

#[test]
fn an_evolving_pipe_refuses_what_it_does_not_build_by_name() {
    let mut model = Model::new();
    let wire = square_wire(&mut model, 1.0);
    let profile = xy_face(&mut model, &[wire]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, 3.0));
    let law = PipeLaw::RotationMinimizing;
    let err = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &law,
        evolve(-1.0, 0.0),
        SKIN,
        T,
    )
    .unwrap_err();
    assert!(err.to_string().contains("naught or more"), "{err}");
    let err = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &law,
        evolve(0.0, 1.0),
        SKIN,
        T,
    )
    .unwrap_err();
    assert!(err.to_string().contains("shrinks to a point"), "{err}");
    let err = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::Fixed,
        evolve(1.0, 1.0),
        SKIN,
        T,
    )
    .unwrap_err();
    assert!(err.to_string().contains("fixed section"), "{err}");
}
