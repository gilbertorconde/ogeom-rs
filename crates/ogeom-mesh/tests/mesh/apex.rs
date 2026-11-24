//! A face that runs to a point where its surface collapses (a cone's apex,
//! a sphere's pole) meshes with the apex as one vertex, shared by the fan
//! round it, and no triangle with two corners at one point in space.
//!
//! Each body is closed by one circle: a disc below and a cone or a dome
//! above, the curved face bounded by the circle and the apex or pole as a
//! degenerate edge, so its chart closes against the row where the surface
//! collapses. The half cone runs to the
//! apex along two rulings, the apex a vertex of its outline.
//!
//! Each triangle's corners hold their own parameters, so the mean of a
//! triangle's three lifted onto the surface is a point of the triangle's
//! own patch, and each corner's normal faces the way its triangle does,
//! the tip's included.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::{FastMap, Tolerances};
use ogeom_geom::{CircleCurve, ConeSurface, Curve, LineCurve, PlaneSurface, SphereSurface};
use ogeom_geom::{Surface as _, SurfaceGeometry};
use ogeom_math::{Circle, Cone, Direction, Frame, Plane, Point, Sphere, Vector};
use ogeom_mesh::{Deflection, face_meshes_for};
use ogeom_topo::{EdgeData, Model, Shape, Triangulation, VertexData};

const T: Tolerances = Tolerances::millimetres();

/// The chords the faces are meshed at: fine, middling and coarse for a part
/// a few tens of millimetres across.
const CHORDS: [f64; 3] = [0.01 / 3.0, 0.1, 1.0];

/// A solid of a disc of `radius` in the XY plane and a cap above it, the
/// two sharing the circle that bounds both. `cap` builds the cap's face
/// over that circle, as a file's reader does for a face that names the
/// circle alone.
fn closed_by_a_circle(
    model: &mut Model,
    radius: f64,
    cap: impl FnOnce(&mut Model, &Shape) -> Shape,
) -> Shape {
    let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, radius, T).unwrap()).into();
    let rim = ogeom_algo::make_edge(model, circle, (0.0, std::f64::consts::TAU), T)
        .unwrap()
        .shape;
    let down = Frame::new(Point::ORIGIN, -Direction::Z, Direction::X, T).unwrap();
    let disc = ogeom_algo::make_face_with_pcurves(
        model,
        PlaneSurface::new(Plane::new(down)).into(),
        &[vec![rim.reversed()]],
        T,
    )
    .unwrap()
    .shape;
    let cap = cap(model, &rim);
    let shell = ogeom_algo::make_shell(model, &[disc, cap]).unwrap().shape;
    ogeom_algo::make_solid(model, &[shell]).unwrap().shape
}

/// A cone of base radius 10 rising to its apex 30 above, its axis pointing
/// down from the apex as the surface is stated, over the wide height range
/// a file's unbounded surface is read with. The apex is a degenerate edge,
/// and the cone's chart closes against the row where it collapses.
fn cone_to_a_point(model: &mut Model) -> (Shape, Point) {
    let down = Frame::new(Point::ORIGIN, -Direction::Z, Direction::X, T).unwrap();
    let cone = Cone::new(down, 10.0, (1.0_f64 / 3.0).atan(), T).unwrap();
    let surface: SurfaceGeometry = ConeSurface::new(cone, (-1e5, 1e5)).unwrap().into();
    let body = closed_by_a_circle(model, 10.0, |model, rim| {
        ogeom_algo::make_apex_band(model, &surface, rim, T).unwrap()
    });
    (body, Point::new(0.0, 0.0, 30.0))
}

/// A hemisphere of radius 20 over its equator, its pole at the top a
/// degenerate edge.
fn dome(model: &mut Model) -> (Shape, Point) {
    let pole = Point::new(0.0, 0.0, 20.0);
    let surface: SurfaceGeometry =
        SphereSurface::new(Sphere::new(Frame::WORLD, 20.0, T).unwrap()).into();
    let body = closed_by_a_circle(model, 20.0, |model, rim| {
        let vertex = model.add_vertex(VertexData::new(pole));
        let mut data = EdgeData::new();
        data.degenerate = true;
        let collapsed = model.add_edge(data, &[vertex.clone(), vertex]).unwrap();
        ogeom_algo::make_revolution_band(model, &surface, rim, &collapsed, T).unwrap()
    });
    (body, pole)
}

/// Half a cone of base radius 10 and height 30, cut along the XZ plane:
/// bounded by half the base circle and the two rulings from its ends to
/// the apex, and closed by a half disc and the triangle across the cut.
fn half_cone(model: &mut Model) -> (Shape, Point) {
    let apex = Point::new(0.0, 0.0, 30.0);
    let (east, west) = (Point::new(10.0, 0.0, 0.0), Point::new(-10.0, 0.0, 0.0));
    let [top, e, w] = [apex, east, west].map(|p| ogeom_algo::make_vertex(model, p).shape);
    let mut between = |curve: Curve, range: (f64, f64), from: &Shape, to: &Shape| {
        ogeom_algo::make_edge_between(model, curve, range, from, to, T)
            .unwrap()
            .shape
    };
    let line = |a: Point, b: Point| -> Curve { LineCurve::segment(a, b, T).unwrap().into() };
    let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 10.0, T).unwrap()).into();
    let rim = between(circle, (0.0, std::f64::consts::PI), &e, &w);
    let up_e = between(line(east, apex), (0.0, east.distance(apex)), &e, &top);
    let up_w = between(line(west, apex), (0.0, west.distance(apex)), &w, &top);
    let across = between(line(west, east), (0.0, 20.0), &w, &e);

    let mut face = |surface: SurfaceGeometry, edges: Vec<Shape>| {
        ogeom_algo::make_face_with_pcurves(model, surface, &[edges], T)
            .unwrap()
            .shape
    };
    let down = Frame::new(apex, -Direction::Z, Direction::X, T).unwrap();
    let cone = Cone::new(down, 0.0, (1.0_f64 / 3.0).atan(), T).unwrap();
    let side = face(
        ConeSurface::new(cone, (0.0, 1e3)).unwrap().into(),
        vec![rim.clone(), up_w.clone(), up_e.reversed()],
    );
    let below = Frame::new(Point::ORIGIN, -Direction::Z, Direction::X, T).unwrap();
    let floor = face(
        PlaneSurface::new(Plane::new(below)).into(),
        vec![rim.reversed(), across.reversed()],
    );
    let cut = Frame::new(Point::ORIGIN, -Direction::Y, Direction::X, T).unwrap();
    let wall = face(
        PlaneSurface::new(Plane::new(cut)).into(),
        vec![across, up_e, up_w.reversed()],
    );
    let shell = ogeom_algo::make_shell(model, &[side, floor, wall])
        .unwrap()
        .shape;
    (ogeom_algo::make_solid(model, &[shell]).unwrap().shape, apex)
}

/// The triangles with two corners within the confusion distance of each
/// other, by their corners.
fn collapsed_triangles(mesh: &Triangulation) -> Vec<[Point; 3]> {
    mesh.triangles
        .iter()
        .map(|t| t.map(|i| mesh.positions[i as usize]))
        .filter(|[a, b, c]| {
            a.distance(*b) <= T.confusion()
                || b.distance(*c) <= T.confusion()
                || c.distance(*a) <= T.confusion()
        })
        .collect()
}

/// The distance from `p` to the triangle `abc`.
fn to_triangle(p: Point, [a, b, c]: [Point; 3]) -> f64 {
    let n = (b - a).cross(c - a);
    let area2 = n.dot(n);
    if area2 > 0.0 {
        let inside = [(a, b), (b, c), (c, a)]
            .iter()
            .all(|&(x, y)| (y - x).cross(p - x).dot(n) >= 0.0);
        if inside {
            return (p - a).dot(n).abs() / area2.sqrt();
        }
    }
    [(a, b), (b, c), (c, a)]
        .iter()
        .map(|&(x, y)| {
            let d = y - x;
            let l2 = d.dot(d);
            let t = if l2 > 0.0 {
                ((p - x).dot(d) / l2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            (x + d * t).distance(p)
        })
        .fold(f64::INFINITY, f64::min)
}

/// The triangles of `mesh` on `surface` whose corners' mean parameters lift
/// to a point off the triangle's patch, or whose corner at `apex` has a
/// normal facing away from the triangle: the triangle's index and why.
///
/// On its patch is within the chord of the triangle's plane and within
/// the triangle's longest side of the triangle: a corner at the apex
/// holding a `u` from the other side of the seam puts the mean on the far
/// side of the cone, off the plane by a good part of its radius there.
fn off_their_patch(
    mesh: &Triangulation,
    surface: &SurfaceGeometry,
    chord: f64,
    apex: Point,
) -> Vec<(usize, String)> {
    let mut bad = Vec::new();
    for (i, t) in mesh.triangles.iter().enumerate() {
        let corners = t.map(|k| mesh.positions[k as usize]);
        let [a, b, c] = corners;
        let uv = mesh
            .triangle_parameters(i)
            .expect("a face's mesh holds its corners' parameters");
        let mean = (
            (uv[0].0 + uv[1].0 + uv[2].0) / 3.0,
            (uv[0].1 + uv[1].1 + uv[2].1) / 3.0,
        );
        let on = surface.point_at(mean.0, mean.1, T).unwrap();
        let facing: Vector = (b - a).cross(c - a);
        let Ok(normal) = facing.normalized(T) else {
            continue;
        };
        let off = (on - a).dot(normal).abs();
        let longest = a.distance(b).max(b.distance(c)).max(c.distance(a));
        if off > chord + T.confusion() || to_triangle(on, corners) > longest {
            bad.push((i, format!("the mean of {uv:?} lifts {off} off its plane")));
        }
        let normals = mesh
            .triangle_normals(i)
            .expect("a face's mesh holds its corners' normals");
        for k in 0..3 {
            if corners[k].distance(apex) <= T.confusion() && normals[k].dot(facing) <= 0.0 {
                bad.push((i, format!("the tip's normal {:?} faces away", normals[k])));
            }
        }
    }
    bad
}

/// Whether the meshes, welded by position, close: every side walked once
/// each way.
fn closes(meshes: &[&Triangulation]) -> Result<(), String> {
    let key = |p: Point| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
    let mut walked: FastMap<([u64; 3], [u64; 3]), i32> = FastMap::default();
    for mesh in meshes {
        for t in &mesh.triangles {
            for k in 0..3 {
                let (a, b) = (
                    key(mesh.positions[t[k] as usize]),
                    key(mesh.positions[t[(k + 1) % 3] as usize]),
                );
                if a < b {
                    *walked.entry((a, b)).or_default() += 1;
                } else {
                    *walked.entry((b, a)).or_default() -= 1;
                }
            }
        }
    }
    let open: Vec<_> = walked.iter().filter(|(_, n)| **n != 0).collect();
    if open.is_empty() {
        Ok(())
    } else {
        Err(format!("{} sides not walked once each way", open.len()))
    }
}

/// Mesh `build`'s body at each of [`CHORDS`] and hold every face to: no
/// triangle with two corners at one point, at most one vertex at the apex,
/// the chord met, and the faces welded by position closing.
fn check(name: &str, build: fn(&mut Model) -> (Shape, Point)) {
    let mut model = Model::new();
    let (body, apex) = build(&mut model);
    for chord in CHORDS {
        let deflection = Deflection {
            chord,
            ..Deflection::default()
        };
        let (faces, _) = face_meshes_for(&model, &body, deflection, None, T).unwrap();
        for (i, (face, mesh)) in faces.iter().enumerate() {
            let data = model.node(face).unwrap().data().as_face().unwrap();
            let surface = model.geometry().surface(data.surface).unwrap();
            let bad = off_their_patch(mesh.as_ref().unwrap(), surface, chord, apex);
            assert!(
                bad.is_empty(),
                "{name} at chord {chord}: face {i} has {} triangles off their patch: {:?}",
                bad.len(),
                &bad[..bad.len().min(4)]
            );
        }
        let meshes: Vec<&Triangulation> = faces
            .iter()
            .map(|(_, mesh)| mesh.as_ref().unwrap())
            .collect();
        let mut reached = false;
        for (i, mesh) in meshes.iter().enumerate() {
            let collapsed = collapsed_triangles(mesh);
            assert!(
                collapsed.is_empty(),
                "{name} at chord {chord}: face {i} has {} triangles with two corners at one point: {collapsed:?}",
                collapsed.len()
            );
            let at_apex = mesh
                .positions
                .iter()
                .filter(|p| p.distance(apex) <= T.confusion())
                .count();
            assert!(
                at_apex <= 1,
                "{name} at chord {chord}: face {i} has {at_apex} vertices at the apex"
            );
            reached |= at_apex == 1;
            assert!(
                mesh.deflection_met,
                "{name} at chord {chord}: face {i} misses the chord"
            );
        }
        assert!(reached, "{name} at chord {chord}: no face reaches the apex");
        if let Err(why) = closes(&meshes) {
            panic!("{name} at chord {chord}: {why}");
        }
    }
}

#[test]
fn a_cone_closing_on_its_apex_meshes_the_apex_as_one_vertex() {
    check("cone", cone_to_a_point);
}

#[test]
fn a_dome_closing_on_its_pole_meshes_the_pole_as_one_vertex() {
    check("dome", dome);
}

#[test]
fn a_half_cone_bounded_at_its_apex_meshes_it_as_one_vertex() {
    check("half cone", half_cone);
}
