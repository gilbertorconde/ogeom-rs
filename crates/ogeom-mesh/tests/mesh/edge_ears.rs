//! Two quarter cylinders of one radius at right angles, mitred where they
//! meet: a block's two rounded edges running into one corner. They share a
//! quarter ellipse, and at a coarse chord its polyline bends out of the
//! chart's straight chords in both cylinders' charts. Neither face may draw
//! a triangle across three of the ellipse's points: the face across could
//! draw the same one, and welded, the chord under the two would be four
//! triangles' side.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::{FastMap, Tolerances};
use ogeom_geom::{CircleCurve, CylinderSurface, EllipseCurve, LineCurve};
use ogeom_math::{Axis, Circle, Cylinder, Direction, Ellipse, Frame, Point, Vector};
use ogeom_mesh::{Deflection, face_meshes_for};
use ogeom_topo::{Model, Shape, Triangulation};

const T: Tolerances = Tolerances::millimetres();

fn direction(x: f64, y: f64, z: f64) -> Direction {
    Direction::new(Vector::new(x, y, z), T).unwrap()
}

fn frame(origin: Point, z: Direction, x: Direction) -> Frame {
    Frame::new(origin, z, x, T).unwrap()
}

fn line(model: &mut Model, from: &Shape, to: &Shape, a: Point, b: Point) -> Shape {
    let along = b - a;
    let axis = Axis::new(a, Direction::new(along, T).unwrap());
    let curve = LineCurve::new(axis).into();
    ogeom_algo::make_edge_between(model, curve, (0.0, along.magnitude()), from, to, T)
        .unwrap()
        .shape
}

/// A shell of the two quarter cylinders, radius 4: one along `y` about
/// `x = 8, z = 14`, from `y = -8` to the mitre, and one along `x` about
/// `y = 9, z = 14`, from `x = -20` to the mitre. The mitre is the plane
/// `x - y = -1`, through the corner where the block's top edges met.
fn mitred() -> (Model, Shape) {
    let mut model = Model::new();
    let m = &mut model;
    let vertex = |m: &mut Model, x, y, z| ogeom_algo::make_vertex(m, Point::new(x, y, z)).shape;
    let top = vertex(m, 12.0, 13.0, 14.0);
    let bottom = vertex(m, 8.0, 9.0, 10.0);
    let (near_top, near_bottom) = (vertex(m, 12.0, -8.0, 14.0), vertex(m, 8.0, -8.0, 10.0));
    let (far_top, far_bottom) = (vertex(m, -20.0, 13.0, 14.0), vertex(m, -20.0, 9.0, 10.0));

    // The mitre's quarter ellipse, from the top corner down to the bottom
    // of both rounds.
    let half = std::f64::consts::FRAC_1_SQRT_2;
    let ellipse = Ellipse::new(
        frame(
            Point::new(8.0, 9.0, 14.0),
            direction(-half, half, 0.0),
            direction(-half, -half, 0.0),
        ),
        4.0 * std::f64::consts::SQRT_2,
        4.0,
        T,
    )
    .unwrap();
    let pi = std::f64::consts::PI;
    let mitre = ogeom_algo::make_edge_between(
        m,
        EllipseCurve::new(ellipse).into(),
        (pi, 1.5 * pi),
        &top,
        &bottom,
        T,
    )
    .unwrap()
    .shape;

    // Along y.
    let crest = line(
        m,
        &near_top,
        &top,
        Point::new(12.0, -8.0, 14.0),
        Point::new(12.0, 13.0, 14.0),
    );
    let foot = line(
        m,
        &bottom,
        &near_bottom,
        Point::new(8.0, 9.0, 10.0),
        Point::new(8.0, -8.0, 10.0),
    );
    let near_end = Circle::new(
        frame(Point::new(8.0, -8.0, 14.0), Direction::Y, Direction::X),
        4.0,
        T,
    )
    .unwrap();
    let near_end = ogeom_algo::make_edge_between(
        m,
        CircleCurve::new(near_end).into(),
        (0.0, 0.5 * pi),
        &near_top,
        &near_bottom,
        T,
    )
    .unwrap()
    .shape;
    let along_y = CylinderSurface::new(
        Cylinder::new(
            frame(
                Point::new(8.0, 0.0, 14.0),
                Direction::Y,
                direction(half, 0.0, -half),
            ),
            4.0,
            T,
        )
        .unwrap(),
        (-100.0, 100.0),
    )
    .unwrap();
    let first = ogeom_algo::make_face_with_pcurves(
        m,
        along_y.into(),
        &[vec![mitre.clone(), foot, near_end.reversed(), crest]],
        T,
    )
    .unwrap()
    .shape
    .reversed();

    // Along x.
    let crest = line(
        m,
        &top,
        &far_top,
        Point::new(12.0, 13.0, 14.0),
        Point::new(-20.0, 13.0, 14.0),
    );
    let foot = line(
        m,
        &far_bottom,
        &bottom,
        Point::new(-20.0, 9.0, 10.0),
        Point::new(8.0, 9.0, 10.0),
    );
    let far_end = Circle::new(
        frame(Point::new(-20.0, 9.0, 14.0), -Direction::X, Direction::Y),
        4.0,
        T,
    )
    .unwrap();
    let far_end = ogeom_algo::make_edge_between(
        m,
        CircleCurve::new(far_end).into(),
        (0.0, 0.5 * pi),
        &far_top,
        &far_bottom,
        T,
    )
    .unwrap()
    .shape;
    let along_x = CylinderSurface::new(
        Cylinder::new(
            frame(
                Point::new(2.5, 9.0, 14.0),
                -Direction::X,
                direction(0.0, half, -half),
            ),
            4.0,
            T,
        )
        .unwrap(),
        (-100.0, 100.0),
    )
    .unwrap();
    let second = ogeom_algo::make_face_with_pcurves(
        m,
        along_x.into(),
        &[vec![foot, mitre.reversed(), crest, far_end]],
        T,
    )
    .unwrap()
    .shape
    .reversed();

    let shell = ogeom_algo::make_shell(m, &[first, second]).unwrap().shape;
    (model, shell)
}

/// How many triangles use each side of `mesh`, by position.
fn sides(mesh: &Triangulation) -> FastMap<(u32, u32), usize> {
    let mut uses: FastMap<(u32, u32), usize> = FastMap::default();
    for t in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            *uses.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    uses
}

#[test]
fn a_mitre_meshes_no_triangle_across_its_ellipse_in_either_face() {
    let (model, shell) = mitred();
    for chord in [0.0125, 0.05, 0.1, 0.25, 0.5] {
        let deflection = Deflection {
            chord,
            ..Deflection::default()
        };
        let (faces, _) = face_meshes_for(&model, &shell, deflection, None, T).unwrap();
        assert_eq!(faces.len(), 2);
        let mut joined = Triangulation::new();
        for (n, (_, mesh)) in faces.iter().enumerate() {
            let mesh = mesh.as_ref().unwrap();
            assert!(mesh.deflection_met, "face {n} at {chord}: misses its chord");
            let uses = sides(mesh);
            let border = |a: u32, b: u32| uses.get(&(a.min(b), a.max(b))) == Some(&1);
            let on_mitre = |i: u32| {
                let p = mesh.positions[i as usize];
                (p.x - p.y + 1.0).abs() <= 1e-9
            };
            for t in &mesh.triangles {
                assert!(
                    !t.iter().all(|&i| on_mitre(i))
                        || !(0..3).any(|k| border(t[k], t[(k + 1) % 3])
                            && border(t[(k + 1) % 3], t[(k + 2) % 3])),
                    "face {n} at {chord}: {:?} spans the ellipse",
                    t.map(|i| mesh.positions[i as usize])
                );
            }
            joined.append(mesh);
        }
        let welded = joined.welded(T);
        let more = sides(&welded).into_iter().filter(|&(_, n)| n > 2).count();
        assert_eq!(
            more, 0,
            "at {chord}: {more} sides of more than two triangles"
        );
    }
}
