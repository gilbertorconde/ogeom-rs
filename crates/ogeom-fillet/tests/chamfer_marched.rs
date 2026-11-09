//! Chamfers on edges only the marching machinery seats: a bore's rim on a
//! drum.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{Chamfer, chamfer_edges_with};
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// A drum of radius 10 and height 30 with a bore of radius 3 across it
/// along y at height 15.
fn holed_drum(m: &mut Model) -> Shape {
    let drum = ogeom_algo::make_cylinder(m, at(0.0, 0.0, 0.0), 10.0, 30.0, T)
        .unwrap()
        .shape;
    let f = Frame::new(Point::new(0.0, -20.0, 15.0), Direction::Y, Direction::X, T).unwrap();
    let bore = ogeom_algo::make_cylinder(m, f, 3.0, 40.0, T).unwrap().shape;
    ogeom_bool::cut(m, &drum, &bore, T).unwrap().shape
}

/// The edge of `s` passing nearest `p`.
fn edge_near(m: &Model, s: &Shape, p: [f64; 3]) -> Shape {
    let q = Point::new(p[0], p[1], p[2]);
    let d = |e: &Shape| {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
            return f64::INFINITY;
        };
        (0..=64)
            .map(|i| {
                c.point_at(a + (b - a) * f64::from(i) / 64.0, T)
                    .unwrap()
                    .distance(q)
            })
            .fold(f64::INFINITY, f64::min)
    };
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| d(a).total_cmp(&d(b)))
        .unwrap()
}

/// How far `p` stands off the drum's wall (radius 10 about z) and off the
/// bore's (radius 3 about the line x = 0, z = 15).
fn off_drum(p: Point) -> f64 {
    (p.x.hypot(p.y) - 10.0).abs()
}

fn off_bore(p: Point) -> f64 {
    (p.x.hypot(p.z - 15.0) - 3.0).abs()
}

/// The nearest point of the crease to `p`, refined on the crease's own
/// parameter: where `p` stands in the plane square to the crease; and the
/// crease's unit tangent there.
fn crease_near(p: Point) -> (Point, Vector) {
    let at = |t: f64| {
        let x = 3.0 * t.cos();
        Point::new(x, (100.0 - x * x).sqrt(), 3.0f64.mul_add(t.sin(), 15.0))
    };
    let n = 4096;
    let mut best = (f64::INFINITY, 0.0);
    for i in 0..n {
        let t = core::f64::consts::TAU * f64::from(i) / f64::from(n);
        let d = at(t).distance(p);
        if d < best.0 {
            best = (d, t);
        }
    }
    let step = core::f64::consts::TAU / f64::from(n);
    let (mut lo, mut hi) = (best.1 - step, best.1 + step);
    for _ in 0..80 {
        let (a, b) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
        if at(a).distance(p) < at(b).distance(p) {
            hi = b;
        } else {
            lo = a;
        }
    }
    let t = f64::midpoint(lo, hi);
    let h = 1e-6;
    let tangent = at(t + h) - at(t - h);
    (at(t), tangent / tangent.magnitude())
}

fn nearest_on_crease(p: Point) -> Point {
    crease_near(p).0
}

/// The holed drum's bore rim chamfered by `spec` (handed the drum's face
/// and the bore's), checked: back come the bevel's rulings at 32
/// stations, each from its drum end to its bore end.
fn bore_rulings(spec: impl Fn(&Shape, &Shape) -> Chamfer) -> Vec<(Point, Point)> {
    let mut m = Model::new();
    let s = holed_drum(&mut m);
    let rim = edge_near(&m, &s, [0.0, 9.5, 18.0]);
    let beside: Vec<Shape> = explore_unique(&m, &s, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            explore_unique(&m, f, ShapeType::Edge)
                .unwrap()
                .iter()
                .any(|x| x.is_same(&rim))
        })
        .collect();
    let surface = |m: &Model, f: &Shape| {
        m.node(f)
            .and_then(|n| n.data().as_face())
            .and_then(|d| m.geometry().surface(d.surface))
            .cloned()
            .unwrap()
    };
    let is_drum = |f: &Shape| match surface(&m, f) {
        ogeom_geom::SurfaceGeometry::Cylinder(c) => (c.cylinder().radius() - 10.0).abs() < 1e-9,
        _ => false,
    };
    let drum = beside
        .iter()
        .find(|f| is_drum(f))
        .expect("the drum")
        .clone();
    let bore = beside
        .iter()
        .find(|f| !is_drum(f))
        .expect("the bore")
        .clone();
    let built = chamfer_edges_with(&mut m, &s, &[(rim.clone(), spec(&drum, &bore))], T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(built.history.is_deleted(&rim));
    let bevel = explore_unique(&m, &built.shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| matches!(surface(&m, f), ogeom_geom::SurfaceGeometry::BSpline(_)))
        .expect("the bevel");
    let geometry = surface(&m, &bevel);
    let ((u0, u1), (v0, v1)) = geometry.domain();
    (0..32)
        .map(|i| {
            let v = (v1 - v0).mul_add(f64::from(i) / 32.0, v0);
            let (a, b) = (
                geometry.point_at(u0, v, T).unwrap(),
                geometry.point_at(u1, v, T).unwrap(),
            );
            if off_drum(a) < off_drum(b) {
                (a, b)
            } else {
                (b, a)
            }
        })
        .collect()
}

/// The band's fit target: how far a rail may stand off its face.
const FIT: f64 = 2e-4;

/// A leg of `length` measured along a section of curvature at most `kappa`
/// spans a chord between `length (1 - (kappa length)² / 24)` and `length`.
fn leg_fits(chord: f64, length: f64, kappa: f64) {
    let shortest = length * (1.0 - (kappa * length).powi(2) / 24.0);
    assert!(
        chord <= length + FIT && chord >= shortest - FIT,
        "a chord of {chord} for a leg of {length}"
    );
}

/// How sharply a cylinder of radius `r`'s section by the plane square to
/// the crease bends within `length` of the crease. The crease's tangent
/// lies in the face, so at the crease that plane holds the face's normal: a
/// normal section, of curvature at most `1 / r`. Along the leg the normal
/// turns by at most `length / r`, which leans the plane off it by as much,
/// and by Meusnier the section's curvature is the normal one over the
/// lean's cosine.
fn kappa(r: f64, length: f64) -> f64 {
    1.0 / (r * (length / r).cos())
}

#[test]
fn a_bore_rim_on_a_drum_chamfers_symmetric() {
    for (on_drum, on_bore) in bore_rulings(|_, _| Chamfer::Symmetric(0.5)) {
        assert!(off_drum(on_drum) <= FIT && off_bore(on_bore) <= FIT);
        leg_fits(
            on_drum.distance(nearest_on_crease(on_drum)),
            0.5,
            kappa(10.0, 0.5),
        );
        leg_fits(
            on_bore.distance(nearest_on_crease(on_bore)),
            0.5,
            kappa(3.0, 0.5),
        );
    }
}

#[test]
fn a_bore_rim_on_a_drum_chamfers_by_two_distances() {
    let rulings = bore_rulings(|drum, _| Chamfer::Distances {
        face: drum.clone(),
        on_face: 0.8,
        on_other: 0.4,
    });
    for (on_drum, on_bore) in rulings {
        assert!(off_drum(on_drum) <= FIT && off_bore(on_bore) <= FIT);
        leg_fits(
            on_drum.distance(nearest_on_crease(on_drum)),
            0.8,
            kappa(10.0, 0.8),
        );
        leg_fits(
            on_bore.distance(nearest_on_crease(on_bore)),
            0.4,
            kappa(3.0, 0.4),
        );
    }
}

#[test]
fn a_bore_rim_on_a_drum_chamfers_by_distance_and_angle() {
    let angle = 40.0_f64.to_radians();
    let rulings = bore_rulings(|drum, _| Chamfer::Angle {
        face: drum.clone(),
        distance: 0.6,
        angle,
    });
    for (on_drum, on_bore) in rulings {
        assert!(off_drum(on_drum) <= FIT && off_bore(on_bore) <= FIT);
        leg_fits(
            on_drum.distance(nearest_on_crease(on_drum)),
            0.6,
            kappa(10.0, 0.6),
        );
        // The ruling lies in the plane square to the crease through the
        // drum's contact, and leaves the drum's section there at the angle:
        // its sine against the drum's normal within that plane.
        let (_, along) = crease_near(on_drum);
        let ruling = on_bore - on_drum;
        assert!(ruling.dot(along).abs() <= FIT, "{}", ruling.dot(along));
        let normal = Vector::new(on_drum.x, on_drum.y, 0.0) / 10.0;
        let normal = normal - along * normal.dot(along);
        let normal = normal / normal.magnitude();
        let off = (ruling.dot(normal) / ruling.magnitude()).abs().asin();
        assert!((off - angle).abs() <= 1e-3, "{} degrees", off.to_degrees());
    }
}
