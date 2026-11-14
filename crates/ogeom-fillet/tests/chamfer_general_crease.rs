//! A chamfer on a crease no closed form seats: where a tilted plane meets
//! a cylinder in an ellipse, against the volume its section sweeps along
//! the crease.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::chamfer_edges;
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
/// The cap's tilt: it stands at z = 10 + TILT y.
const TILT: f64 = 0.3;

fn up_at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn vol(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(
        m,
        s,
        ogeom_mesh::Deflection {
            chord: 1e-5,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

/// A half cylinder of radius 3 (x >= 0) from z = 0 up to the tilted cap
/// z = 10 + 0.3 y: the cap meets the wall along half an ellipse, square to
/// the flat x = 0 at both ends.
fn tilted_d_block(m: &mut Model) -> Shape {
    let full = ogeom_algo::make_cylinder(m, up_at(0.0, 0.0, 0.0), 3.0, 14.0, T)
        .unwrap()
        .shape;
    let half = ogeom_algo::make_box(m, up_at(-4.0, -4.0, -1.0), (4.0, 8.0, 16.0), T)
        .unwrap()
        .shape;
    let block = ogeom_bool::cut(m, &full, &half, T).unwrap().shape;
    // Everything above the tilted cap, by a box standing on its plane.
    let n = Vector::new(0.0, -TILT, 1.0);
    let z = n / n.magnitude();
    let x = Vector::X;
    let corner = Point::new(-30.0, 0.0, 10.0) - z.cross(x) * 30.0;
    let frame = Frame::new(corner, Direction::new(z, T).unwrap(), Direction::X, T).unwrap();
    let lid = ogeom_algo::make_box(m, frame, (60.0, 60.0, 30.0), T)
        .unwrap()
        .shape;
    ogeom_bool::cut(m, &block, &lid, T).unwrap().shape
}

/// The crease at `phi`, and its derivative.
fn crease(phi: f64) -> (Vector, Vector) {
    let (s, c) = phi.sin_cos();
    (
        Vector::new(3.0 * c, 3.0 * s, TILT.mul_add(3.0 * s, 10.0)),
        Vector::new(-3.0 * s, 3.0 * c, TILT * 3.0 * c),
    )
}

/// The section of the chamfer `d` at the crease's `phi`, in the plane
/// square to the crease: its area and its centroid's offset from the
/// crease point. The cap's leg runs straight along the cap, square to the
/// crease, into the cap; the wall's along the wall's section, `d` of arc
/// down, walked in fine steps held to the wall and the plane, their chords
/// summed. The region
/// is the triangle of the crease point and the two contacts with the
/// segment between the chord to the wall's contact and the wall.
fn section(phi: f64, d: f64) -> (f64, Vector) {
    let (p, dp) = crease(phi);
    let t = dp / dp.magnitude();
    let cap_n = Vector::new(0.0, -TILT, 1.0) / Vector::new(0.0, -TILT, 1.0).magnitude();
    let mut into_cap = cap_n.cross(t);
    if into_cap.dot(Vector::new(-p.x, -p.y, 0.0)) < 0.0 {
        into_cap = -into_cap;
    }
    let c1 = p + into_cap * d;
    // The wall's section: steps in the plane square to `t`, held to the
    // cylinder, downward.
    let h = d / 4000.0;
    let mut q = p;
    let mut arc = vec![p];
    let mut walked = 0.0;
    while walked < d {
        let radial = Vector::new(q.x, q.y, 0.0) / q.x.hypot(q.y);
        let mut dir = t.cross(radial);
        if dir.z > 0.0 {
            dir = -dir;
        }
        let mut next = q + dir * h;
        for _ in 0..8 {
            next = next - t * (next - p).dot(t);
            let r = next.x.hypot(next.y);
            next = Vector::new(next.x * 3.0 / r, next.y * 3.0 / r, next.z);
        }
        let step = (next - q).magnitude();
        if walked + step >= d {
            // The last step cut back to the leg's length along its chord.
            next = q + (next - q) * ((d - walked) / step);
            walked = d;
        } else {
            walked += step;
        }
        q = next;
        arc.push(q);
    }
    // The polygon: the crease point, the cap's contact, the wall's, then
    // back up the wall's arc.
    let mut ring = vec![p, c1];
    ring.extend(arc.iter().rev());
    let e1 = into_cap;
    let e2 = t.cross(e1);
    let flat: Vec<(f64, f64)> = ring
        .iter()
        .map(|x| ((*x - p).dot(e1), (*x - p).dot(e2)))
        .collect();
    let (mut a, mut cx, mut cy) = (0.0, 0.0, 0.0);
    for k in 0..flat.len() {
        let (x0, y0) = flat[k];
        let (x1, y1) = flat[(k + 1) % flat.len()];
        let cross = x0 * y1 - x1 * y0;
        a += cross / 2.0;
        cx += (x0 + x1) * cross / 6.0;
        cy += (y0 + y1) * cross / 6.0;
    }
    (a.abs(), (e1 * cx + e2 * cy) / a)
}

/// The volume the section sweeps along the crease from `phi = -π/2` to
/// `π/2`: `∫ A (1 - κ N . c) ds`, the crease's curvature vector read from
/// its unit tangent's turn, by Simpson's rule.
fn swept(d: f64) -> f64 {
    let n = 400;
    let (lo, hi) = (-core::f64::consts::FRAC_PI_2, core::f64::consts::FRAC_PI_2);
    let step = (hi - lo) / f64::from(n);
    let unit = |phi: f64| {
        let (_, dp) = crease(phi);
        dp / dp.magnitude()
    };
    let f = |phi: f64| {
        let (_, dp) = crease(phi);
        let speed = dp.magnitude();
        let h = 1e-5;
        let curvature = (unit(phi + h) - unit(phi - h)) / (2.0 * h * speed);
        let (area, centroid) = section(phi, d);
        area * (1.0 - curvature.dot(centroid)) * speed
    };
    let mut sum = f(lo) + f(hi);
    for k in 1..n {
        let phi = step.mul_add(f64::from(k), lo);
        sum += f(phi) * if k % 2 == 1 { 4.0 } else { 2.0 };
    }
    sum * step / 3.0
}

/// The tilted cap's edge on the wall chamfered 0.3: valid, the edge
/// deleted into its bevel, a B-spline, its ends flush on the flat; the
/// volume taken off the section's sweep along the crease.
#[test]
fn a_tilted_cap_s_edge_on_a_cylinder_chamfers_to_its_swept_section() {
    let d = 0.3;
    let mut m = Model::new();
    let s = tilted_d_block(&mut m);
    let crease_edge = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (c, _) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
            !matches!(c, ogeom_geom::Curve::Line(_))
                && ogeom_algo::shape_bounds(&m, e, T).unwrap().low().unwrap().z > 8.0
        })
        .expect("the cap's edge on the wall");
    let before = vol(&m, &s);
    let built = chamfer_edges(&mut m, &s, std::slice::from_ref(&crease_edge), d, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    assert!(built.history.is_deleted(&crease_edge));
    assert_eq!(built.history.generated(&crease_edge).len(), 1, "the bevel");
    let lost = before - vol(&m, &built.shape);
    let want = swept(d);
    assert!(
        (lost - want).abs() <= 1e-6 * want,
        "lost {lost} against {want}"
    );
}
