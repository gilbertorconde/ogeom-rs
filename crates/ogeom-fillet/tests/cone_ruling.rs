//! A cone meeting a plane along a ruling, rounded through the engine: the
//! ball's centre runs on the conic the plane's and the cone's offsets cut,
//! with no closed form, so the strip is fitted through its exact sections.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// Each face of `shape`: its placed surface and its trim.
fn faces_of(model: &Model, shape: &Shape) -> Vec<(ogeom_geom::SurfaceGeometry, Shape)> {
    use ogeom_geom::Transformable as _;
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .map(|f| {
            let data = model.node(&f).unwrap().data().as_face().unwrap();
            let surface = model
                .geometry()
                .surface(data.surface)
                .unwrap()
                .clone()
                .transformed(&f.transform(model.datums()).unwrap(), T)
                .unwrap();
            (surface, f)
        })
        .collect()
}

/// How far `p` stands from the nearest of `faces` whose trim holds its foot.
fn off(
    model: &Model,
    faces: &[(ogeom_geom::SurfaceGeometry, Shape)],
    trims: &[ogeom_algo::classify::FaceTrim<'_>],
    p: Point,
) -> f64 {
    let mut best = f64::INFINITY;
    for ((surface, _), trim) in faces.iter().zip(trims) {
        let foot = ogeom_algo::project_on_surface(surface, p, 32, T).unwrap();
        if foot.distance < best
            && trim
                .holds(model, foot.parameters, T.confusion(), T)
                .unwrap()
                .is_inside_or_on()
        {
            best = foot.distance;
        }
    }
    best
}

/// A frustum of radius 6 at its foot and 3 at its top, 6 tall, cut in half
/// through its axis (`y <= 0` kept): the cut face meets the cone along two
/// rulings, each rounded at 1/2. Built through the engine with no boolean,
/// valid, the same solid either way round, and every ball's arc lying on
/// the result within 1e-5: the ball's centre a radius below the cut plane
/// (`y = -1/2`) and a radius in from the cone along its normal, which
/// leans up by the half-angle `a` (`tan a = 1/2`): `|(x, y)| = 6 - z/2 -
/// r / cos a` at height `z`.
#[test]
fn a_half_cone_s_rulings_round_on_their_exact_sections() {
    let r = 0.5;
    let cos_a = 2.0 / 5f64.sqrt();
    let sin_a = 1.0 / 5f64.sqrt();
    let mut volumes = Vec::new();
    for flip in [false, true] {
        let mut model = Model::new();
        let cone = ogeom_algo::make_cone(&mut model, Frame::WORLD, 6.0, 3.0, 6.0, T)
            .unwrap()
            .shape;
        let at = Frame::new(Point::new(-10.0, 0.0, -1.0), Direction::Z, Direction::X, T).unwrap();
        let tool = ogeom_algo::make_box(&mut model, at, (20.0, 10.0, 8.0), T)
            .unwrap()
            .shape;
        let half = ogeom_bool::cut(&mut model, &cone, &tool, T).unwrap().shape;
        let mut rulings: Vec<Shape> = explore_unique(&model, &half, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                ogeom_algo::edge_vertices(&model, e)
                    .unwrap()
                    .is_some_and(|(a, b)| {
                        let p =
                            |v: &Shape| model.node(v).unwrap().data().as_vertex().unwrap().point;
                        let (a, b) = (p(&a), p(&b));
                        a.y.abs() < 1e-9 && b.y.abs() < 1e-9 && (a.z - b.z).abs() > 1.0
                    })
            })
            .collect();
        assert_eq!(rulings.len(), 2);
        if flip {
            rulings.reverse();
        }
        let (built, work) = ogeom_bool::work::counted(|| {
            ogeom_fillet::fillet_edges(&mut model, &half, &rulings, r, T).unwrap()
        });
        assert_eq!(work.face_pairs, 0, "built with no boolean");
        let check = ogeom_algo::check(&model, &built.shape, T).unwrap();
        assert!(check.is_valid(), "{:?}", check.problems);
        volumes.push(
            ogeom_algo::volume_properties(
                &model,
                &built.shape,
                ogeom_mesh::Deflection::default(),
                T,
            )
            .unwrap()
            .mass,
        );
        let faces = faces_of(&model, &built.shape);
        let fine = ogeom_mesh::Deflection {
            chord: T.confusion(),
            ..Default::default()
        };
        let trims: Vec<_> = faces
            .iter()
            .map(|(_, f)| ogeom_algo::classify::FaceTrim::of(&model, f, fine, T).unwrap())
            .collect();
        for side in [1.0_f64, -1.0] {
            for k in 1..20 {
                let z = 6.0 * f64::from(k) / 20.0;
                let rho = 6.0 - 0.5 * z - r / cos_a;
                let x = side * (rho * rho - r * r).sqrt();
                let centre = Point::new(x, -r, z);
                let on_plane = Point::new(x, 0.0, z);
                let radial = Vector::new(x, -r, 0.0) / rho;
                let on_cone = centre + (radial * cos_a + Vector::Z * sin_a) * r;
                let (p, q) = ((on_plane - centre) / r, (on_cone - centre) / r);
                let axis = p.cross(q) / p.cross(q).magnitude();
                let sweep = p.dot(q).clamp(-1.0, 1.0).acos();
                for j in 0..=8 {
                    let t = sweep * f64::from(j) / 8.0;
                    let point = centre + (p * t.cos() + axis.cross(p) * t.sin()) * r;
                    let d = off(&model, &faces, &trims, point);
                    assert!(d < 1e-5, "z {z}: {point:?} is {d} off the result");
                }
            }
        }
    }
    assert!((volumes[0] - volumes[1]).abs() < 1e-9, "either order");
}
