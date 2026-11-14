//! Concave rounds between two walls along their rulings (#972): three
//! cylinders fused side by side, the two concave seams on one side rounded,
//! each strip the exact cylinder whose axis stands a radius off both walls.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// Three cylinders of radius `wall`, `height` tall, their axes `pitch` apart
/// along `x`, fused, and the two concave seams on the side `y < 0`.
fn lobes(model: &mut Model, wall: f64, pitch: f64, height: f64) -> (Shape, Vec<Shape>) {
    let mut lobes = Vec::new();
    for i in 0..3 {
        let frame = Frame::new(
            Point::new(f64::from(i) * pitch, 0.0, 0.0),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        lobes.push(
            ogeom_algo::make_cylinder(model, frame, wall, height, T)
                .unwrap()
                .shape,
        );
    }
    let two = ogeom_bool::fuse(model, &lobes[0], &lobes[1], T)
        .unwrap()
        .shape;
    let solid = ogeom_bool::fuse(model, &two, &lobes[2], T).unwrap().shape;
    let y = -(wall * wall - pitch * pitch / 4.0).sqrt();
    let point = |model: &Model, v: &Shape| {
        v.transform(model.datums())
            .unwrap()
            .apply(model.node(v).unwrap().data().as_vertex().unwrap().point)
    };
    let seams = [0.5 * pitch, 1.5 * pitch]
        .iter()
        .map(|&x| {
            explore_unique(model, &solid, ShapeType::Edge)
                .unwrap()
                .into_iter()
                .find(|e| {
                    ogeom_algo::edge_vertices(model, e)
                        .unwrap()
                        .is_some_and(|(a, b)| {
                            let (a, b) = (point(model, &a), point(model, &b));
                            [a, b]
                                .iter()
                                .all(|p| (p.x - x).abs() < 1e-7 && (p.y - y).abs() < 1e-7)
                                && (a.z - b.z).abs() > 0.5 * height
                        })
                })
                .expect("a concave seam")
        })
        .collect();
    (solid, seams)
}

/// The rounded section's area: three discs, less the two lenses where
/// neighbours overlap, plus each round's fill, the region between the two
/// walls and the ball's arc. A fill is twice its half from the seam's
/// line `x = pitch / 2` out to the contact, under the chord of the wall
/// and above the ball's circle, integrated in closed form.
fn section_area(wall: f64, pitch: f64, radius: f64) -> f64 {
    // The area under `sqrt(r^2 - x^2)` from nought to `x`.
    let under = |r: f64, x: f64| 0.5 * (x * (r * r - x * x).sqrt() + r * r * (x / r).asin());
    let half = 0.5 * pitch;
    let reach = wall + radius;
    // The ball's centre below the seam, and the contact's offset from the
    // wall's axis along `x`.
    let depth = (reach * reach - half * half).sqrt();
    let contact = wall * half / reach;
    let lens =
        2.0 * wall * wall * (half / wall).acos() - pitch * (wall * wall - half * half).sqrt();
    // Between `x = contact` and `x = half` (from the wall's axis): below the
    // wall's lower arc and above the ball's upper arc about `(half,
    // -depth)`, each measured down from `y = 0` as the area over it.
    let below_wall = under(wall, half) - under(wall, contact);
    let below_ball =
        depth * (half - contact) - (under(radius, half - contact) - under(radius, 0.0));
    let fill = 2.0 * (below_ball - below_wall);
    3.0 * core::f64::consts::PI * wall * wall - 2.0 * lens + 2.0 * fill
}

/// #972's two cases: walls 3.58 apart 5.5 rounded at 8, and 4.85 apart 8.5
/// at 4. Built by the engine with no boolean, valid, the same solid either
/// way round, its volume the section's area times the height within 1e-6,
/// and each round the cylinder about the ball's exact centre line.
#[test]
fn concave_rounds_between_fused_cylinders_are_exact_cylinders() {
    for (wall, pitch, height, radius) in [(3.58_f64, 5.5_f64, 14.0, 8.0), (4.85, 8.5, 12.0, 4.0)] {
        let mut volumes = Vec::new();
        for flip in [false, true] {
            let mut model = Model::new();
            let (solid, mut seams) = lobes(&mut model, wall, pitch, height);
            if flip {
                seams.reverse();
            }
            let (built, work) = ogeom_bool::work::counted(|| {
                ogeom_fillet::fillet_edges(&mut model, &solid, &seams, radius, T).unwrap()
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
            let depth = ((wall + radius).powi(2) - 0.25 * pitch * pitch).sqrt();
            for x in [0.5 * pitch, 1.5 * pitch] {
                let centre = Point::new(x, -depth, 0.5 * height);
                let found = explore_unique(&model, &built.shape, ShapeType::Face)
                    .unwrap()
                    .into_iter()
                    .filter(|f| {
                        let data = model.node(f).unwrap().data().as_face().unwrap();
                        let Some(ogeom_geom::SurfaceGeometry::Cylinder(c)) =
                            model.geometry().surface(data.surface)
                        else {
                            return false;
                        };
                        let axis = c.cylinder().axis();
                        let placed = f.transform(model.datums()).unwrap();
                        let (o, d) = (
                            placed.apply(axis.location),
                            placed.apply_vector(axis.direction.vector()),
                        );
                        (c.cylinder().radius() - radius).abs() < T.confusion()
                            && (centre - o).cross(d).magnitude() < T.confusion()
                    })
                    .count();
                assert_eq!(found, 1, "the round at x = {x} is its exact cylinder");
            }
        }
        let want = section_area(wall, pitch, radius) * height;
        for v in &volumes {
            assert!((v - want).abs() < 1e-6, "{v} against {want}");
        }
        assert!((volumes[0] - volumes[1]).abs() < 1e-9, "either order");
    }
}
