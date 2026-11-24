#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
use super::assemble::assemble;
use super::*;
use ogeom_math::{Direction, Frame};
use ogeom_topo::{ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn up_at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// A face of a small box at `(x, y, z)`, 0.2 on a side: its top.
fn probe(m: &mut Model, x: f64, y: f64, z: f64) -> Shape {
    let b = ogeom_algo::make_box(m, up_at(x, y, z), (0.2, 0.2, 0.2), T)
        .unwrap()
        .shape;
    explore_unique(m, &b, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let (p, n) = ogeom_algo::face_normal(m, f, T).unwrap();
            n.z > 0.5 && (p.z - (z + 0.2)).abs() < 1e-9
        })
        .unwrap()
}

fn within(m: &Model, face: &Shape, region: &ogeom_algo::HalfSpaces) -> bool {
    ogeom_algo::certify_face_within(m, face, region, 1e-6, T)
        .unwrap()
        .is_ok()
}

/// The material a face bounds, as the certify reads it: on the D-block
/// (a drum of radius 3 halved along its axis) under the flat (`x >= 0`)
/// and inside the wall; on a tube bored to 2, outside the bore. A face
/// inside is proved within, one past the flat, outside the wall or in
/// the bore is found past.
#[test]
fn a_face_s_material_is_its_side_of_its_plane_or_wall() {
    let mut m = Model::new();
    let full = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 10.0, T)
        .unwrap()
        .shape;
    let half = ogeom_algo::make_box(&mut m, up_at(-4.0, -4.0, -1.0), (4.0, 8.0, 12.0), T)
        .unwrap()
        .shape;
    let d = ogeom_bool::cut(&mut m, &full, &half, T).unwrap().shape;
    let curved = read(&m, &d, T).unwrap().unwrap();
    let of = |kind: fn(&SurfaceGeometry) -> bool| {
        let fi = (0..curved.faces.len())
            .find(|fi| {
                kind(&curved.faces[*fi].surface) && {
                    // The flat at x = 0, not a cap.
                    match &curved.faces[*fi].surface {
                        SurfaceGeometry::Plane(p) => p.plane().normal().vector().x.abs() > 0.5,
                        _ => true,
                    }
                }
            })
            .unwrap();
        material(&curved, fi).unwrap()
    };
    let flat = of(|s| matches!(s, SurfaceGeometry::Plane(_)));
    let wall = of(|s| matches!(s, SurfaceGeometry::Cylinder(_)));
    let inside = probe(&mut m, 1.0, 0.0, 5.0);
    let past_flat = probe(&mut m, -1.0, 0.0, 5.0);
    let past_wall = probe(&mut m, 3.5, 0.0, 5.0);
    assert!(within(&m, &inside, &flat) && within(&m, &inside, &wall));
    assert!(!within(&m, &past_flat, &flat));
    assert!(!within(&m, &past_wall, &wall));

    let outer = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 4.0, T)
        .unwrap()
        .shape;
    let inner = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, -1.0), 2.0, 6.0, T)
        .unwrap()
        .shape;
    let tube = ogeom_bool::cut(&mut m, &outer, &inner, T).unwrap().shape;
    let curved = read(&m, &tube, T).unwrap().unwrap();
    let bore = (0..curved.faces.len())
        .find(|fi| match &curved.faces[*fi].surface {
            SurfaceGeometry::Cylinder(c) => (c.cylinder().radius() - 2.0).abs() < 1e-9,
            _ => false,
        })
        .unwrap();
    let bore = material(&curved, bore).unwrap();
    let in_wall = probe(&mut m, 2.3, 0.0, 2.0);
    let in_bore = probe(&mut m, 0.5, 0.0, 2.0);
    assert!(within(&m, &in_wall, &bore));
    assert!(!within(&m, &in_bore, &bore));
}

/// The D-shaped block's two straight edges rounded at 0.5 by strips a
/// section law would hand over, built here in closed form: the ball
/// between the flat `x = 0` and the wall of radius 3 has its centre at
/// `(r, ±sqrt((3 - r)^2 - r^2))`, touching the flat level with it and
/// the wall on its radius through it; the strip is the cylinder of the
/// ball's radius about that line, its ends arcs about it. The assembly
/// builds them through the caps' flush ends to a valid solid, each strip
/// proved within its faces' material, of the half disc less the
/// line-and-circle corner twice along the 10: the corner the kite from
/// the corner to the contacts and the centre, less the ball's sector,
/// plus the wall's segment beyond the chord.
#[test]
fn round_strips_assemble_through_flush_ends() {
    use core::f64::consts::PI;
    use ogeom_geom::CylinderSurface;
    use ogeom_math::Cylinder;
    let mut m = Model::new();
    let full = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 10.0, T)
        .unwrap()
        .shape;
    let half = ogeom_algo::make_box(&mut m, up_at(-4.0, -4.0, -1.0), (4.0, 8.0, 12.0), T)
        .unwrap()
        .shape;
    let d = ogeom_bool::cut(&mut m, &full, &half, T).unwrap().shape;
    let curved = read(&m, &d, T).unwrap().unwrap();
    let (big, r) = (3.0_f64, 0.5_f64);
    let reach = ((big - r) * (big - r) - r * r).sqrt();
    let mut strips: Vec<Strip> = Vec::new();
    for side in [1.0_f64, -1.0] {
        let edge = *curved
            .edges
            .iter()
            .find(|(_, e)| {
                let (a, b) = (curved.points[&e.ends.0], curved.points[&e.ends.1]);
                a.x.abs() < 1e-9
                    && b.x.abs() < 1e-9
                    && (a.y - big * side).abs() < 1e-9
                    && (a.z - b.z).abs() > 1.0
            })
            .unwrap()
            .0;
        let faces = faces_of(&curved, edge);
        let (cx, cy) = (r, side * reach);
        let on_wall = (cx * big / (big - r), cy * big / (big - r));
        let line = |x: f64, y: f64| -> Curve {
            LineCurve::segment(Point::new(x, y, -5.0), Point::new(x, y, 15.0), T)
                .unwrap()
                .into()
        };
        let contact = |fi: usize| match &curved.faces[fi].surface {
            SurfaceGeometry::Plane(_) => line(0.0, cy),
            _ => line(on_wall.0, on_wall.1),
        };
        let surface: SurfaceGeometry = CylinderSurface::new(
            Cylinder::new(up_at(cx, cy, -5.0), r, T).unwrap(),
            (0.0, 20.0),
        )
        .unwrap()
        .into();
        strips.push(Strip {
            edge,
            faces: [faces[0], faces[1]],
            contacts: [contact(faces[0]), contact(faces[1])],
            surface,
            section: Section::Arc {
                spine: line(cx, cy),
            },
        });
    }
    let mut contacts: FastMap<(usize, Id), Curve> = FastMap::default();
    for st in &strips {
        for k in 0..2 {
            contacts.insert((st.faces[k], st.edge), st.contacts[k].clone());
        }
    }
    let corners = curved_corners(&curved, &contacts, T).unwrap().unwrap();
    let blended: FastSet<Id> = strips.iter().map(|st| st.edge).collect();
    let built = assemble(&mut m, &curved, &strips, &blended, &corners, T)
        .unwrap()
        .expect("the round strips assemble");
    assert!(ogeom_algo::check(&m, &built.shape, T).unwrap().is_valid());
    // The corner: kite less sector plus the wall's segment.
    let (p0, a, c, b) = (
        (0.0, big),
        (0.0, reach),
        (r, reach),
        (r * big / (big - r), reach * big / (big - r)),
    );
    let quad = 0.5
        * ((p0.0 * a.1 - a.0 * p0.1)
            + (a.0 * c.1 - c.0 * a.1)
            + (c.0 * b.1 - b.0 * c.1)
            + (b.0 * p0.1 - p0.0 * b.1))
            .abs();
    let sector = {
        let (u, w) = ((a.0 - c.0, a.1 - c.1), (b.0 - c.0, b.1 - c.1));
        let angle = ((u.0 * w.0 + u.1 * w.1) / (r * r)).clamp(-1.0, 1.0).acos();
        angle * r * r / 2.0
    };
    let segment = {
        let phi = ((b.0 * p0.0 + b.1 * p0.1) / (big * big))
            .clamp(-1.0, 1.0)
            .acos();
        big * big / 2.0 * (phi - phi.sin())
    };
    let corner = quad - sector + segment;
    let want = (PI * big * big / 2.0 - 2.0 * corner) * 10.0;
    let got = ogeom_algo::volume_properties(
        &m,
        &built.shape,
        ogeom_mesh::Deflection {
            chord: 1e-5,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass;
    assert!((got - want).abs() < 1e-6, "volume {got} against {want}");
}
